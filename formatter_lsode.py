from __future__ import annotations
from typing import Iterable
from functools import singledispatchmethod, cached_property

from formatter_fortran import FortranFormatter
from chain import (
    Function, ConstName, Var, NamedFunction, ExplicitFunction, OpaqueFunction,
    ZERO, extract
)

class LsodeFormatter(FortranFormatter):
    """
    Takes FortranFormatter further to allow LSODE use.

    Introduces the concept of variable vector (y), function matrix (F),
    Jacobian matrix and implements function definitions

    Vectors:
    y   solution    Var             the "y" of all ODE systems          INPUT
    f   function    Function        the "f" of y'=f(y)                  INPUT
    c   constants   ConstName       all ConstNames                      EXTRACTED
    a   'all'       NamedFunction   all NamedFunctions (for memoizing)  EXTRACTED

    Note: in the Fortran program, 'c' and 'a' are contiguous with 'y'
    """

    @property
    def y(self) -> tuple[Var, ...]: return self._y

    @cached_property
    def f(self) -> tuple[Function, ...]:
        """
        This encapsulates non-Var Functions in a NamedFunction 'f{i}'
        """
        l = []
        for i, f in enumerate(self._f):
            l.append(f if isinstance(f, Var) else
                        ExplicitFunction(f'f{1+i}', f))
        return tuple(l)

    @property
    def neq(self) -> int: return len(self.y)

    # The below (n_yc and n_yca) are seldom useful, but they're here mostly to
    # make it clear that there is no 'n of y' alone (i.e it is neq which is unambiguous)

    @property
    def n_yc(self) -> int: return self.neq + len(self.c)

    @property
    def n_yca(self) -> int: return self.n_yc + len(self.a_set)

    def __init__(self,
            y: Iterable[Var],
            f: Iterable[Function],
            ) -> None:
        super().__init__()
        self._y = tuple(y)
        self._f = tuple(f)
        if len(self._y) != len(self._f):
            raise ValueError("y and f must have the same length (number of equations)")

    @property
    def jacobian(self) -> tuple[tuple[Function, ...], ...]:
        return tuple([
                tuple([f.der(v) for v in self.y])
                    for f in self.f])

    @cached_property
    def a_set(self) -> set[NamedFunction]:
        """
        Visits everyone in function vector + jacobian to find all NamedFunctions used.

        The result is sorted so that opaque functions come first, for ease of
        finding them to be manually defined later.
        """
        funcs = set(self.f)
        funcs.update({f for row in self.jacobian for f in row})
        return {nf for f in funcs for nf in extract(f, NamedFunction)}

    @property
    def a_human_order(self) -> tuple[NamedFunction, ...]:
        """
        NamedFunction objects in alphabetical order (but opaque and function
        vector first), for ease of manually defining.
        """
        return tuple(sorted(self.a_set,
            key=lambda nf: (1 if isinstance(nf, OpaqueFunction) \
                            else 2 if nf in self.f \
                            else 3,
                            nf.name, len(nf.ders), nf.ders)))

    @cached_property
    def a_eval_order(self) -> tuple[NamedFunction, ...]:
        """
        NamedFunction objects in the correct order of evaluation.

        NamedFunction in the optimized order of evaluation, which allows
        saving each function in a vector to be used by the next functions. That
        is, each is the building block of the next. It is done by counting how
        many functions each uses, and ordering ascending.
        """
        return tuple(sorted(self.a_set,
            key=lambda nf: (len(extract(nf, NamedFunction)),
                            nf.name, len(nf.ders), nf.ders)))

    @cached_property
    def c(self) -> tuple[ConstName, ...]:
        return tuple(sorted(
            {c for f in self.a_set for c in extract(f, ConstName)} ))

    # Repeating needed because @singledispatchmethod doesn't work well with inheritance
    @singledispatchmethod
    def format(self, f: Function) -> str:
        return super().format(f)

    # Overrides for using vectors, e.g turning Var "XG" into "y(1)"

    @format.register
    def _var(self, f: Var) -> str:
        """
        Outputs the Var's reference in the 'y' vector.
        """
        return f'y({1 + self.y.index(f)})'

    @format.register
    def _const_name(self, f: ConstName) -> str:
        """
        Outputs the ConstName's reference in the 'y' vector.

        In LSODE, additional values can be passed in the 'y' vector. We use
        this to pass 'c' and 'a' instead of using global variables.
        *  Fortran 'y' = Python 'y' + 'c' + 'a'
        """
        return f'y({1 + self.neq + self.c.index(f)})'

    @format.register
    def _named(self, f: NamedFunction) -> str:
        """
        Outputs the NamedFunction's reference in the 'y' vector.

        In LSODE, additional values can be passed in the 'y' vector. We use
        this to pass 'c' and 'a' instead of using global variables.
        *  Fortran 'y' = Python 'y' + 'c' + 'a'
        """
        return f'y({1 + self.n_yc + self.a_eval_order.index(f)})'

    def _get_function_name(self, f: NamedFunction) -> str:
        """
        Returns a unique name for Fortran function name.

        Format: <e|o>f_[ders]_<name>
          - ExplicitFunction is 'ef', OpaqueFunction is 'of'
          - [ders] is in the format 'd5_d1' where numbers a Var indices in y
        """
        return '_'.join([
                {ExplicitFunction: 'ef', OpaqueFunction: 'of'}[type(f)],
                f'n{len(extract(f, NamedFunction))}',
                '_'.join([f'd{self.y.index(v)}' for v in f.ders]),
                f.name,
            ])

    def get_function_call(self, f: NamedFunction) -> str:
        return self._get_function_name(f) + '(y)'

    def _get_function_header(self, f: NamedFunction) -> str:
        lines = []
        lines.append(f'double precision function {self.get_function_call(f)} result(res)')
        lines.append(f'  double precision, intent(in) :: y(*)')
        return '\n'.join(lines)

    def _get_entire_function(self, f: NamedFunction, body: str) -> str:
        lines = []
        lines.append(body)
        lines.append('end function')
        return '\n'.join([self._get_function_header(f), '\n'.join(lines)])

    def named_function_definition(self, nf: NamedFunction) -> str:
        body = ('  res = ' + self.format(nf.f)) \
            if isinstance(nf, ExplicitFunction) else '  ! missing opaque definition'
        return self._get_entire_function(nf, body)

    def make_update_a(self):
        lines = []
        lines.append('subroutine update_a(y)')
        lines.append(f'  double precision, intent(inout) :: y(*)')
        lines.append('')
        for nf in self.a_eval_order:
            lines.append(f'  {self.format(nf)} = {self.get_function_call(nf)}')
        lines.append('end subroutine')
        return '\n'.join(lines)

    def make_f(self):
        lines = []
        lines.append('subroutine f(neq, t, y, ydot)')
        lines.append('  integer, intent(in) :: neq')
        lines.append(f'  double precision, intent(in) :: t, y(*), ydot(*)')
        for i, f in enumerate(self.f):
            lines.append(f'  ydot({1+i}) = {self.format(f)}')
        lines.append('end subroutine')
        return '\n'.join(lines)

    def make_jac(self):
        lines = []
        lines.append('subroutine jac(neq, t, ml, mu, pd, nrpd)')
        lines.append('  integer, intent(in) :: neq, ml, mu, nrpd')
        lines.append(f'  double precision, intent(in) :: t, y(*)')
        lines.append(f'  double precision, intent(out) :: pd(nrpd, *)')
        for i, row in enumerate(self.jacobian):
            for j, f in enumerate(row):
                if f != ZERO:
                    lines.append(f'  pd({1+i},{1+j}) = {self.format(f)}')
        lines.append('end subroutine')
        return '\n'.join(lines)
