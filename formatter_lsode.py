from __future__ import annotations
from typing import Iterable
from functools import singledispatchmethod, cached_property

from formatter_fortran import FortranFormatter
from chain import (
    Function, Const, ConstName, Var, NamedFunction, ExplicitFunction,
    OpaqueFunction, ZERO
)
import lsode_program

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
        This encapsulates non-trivial Functions in a NamedFunction 'f{i}'
        """
        l = []
        for i, f in enumerate(self._f):
            l.append(f if isinstance(f, NamedFunction) or _is_trivial(f)
                        else ExplicitFunction(f'f{1+i}', f))
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
                tuple([_unpack_if_trivial(f.der(v)) for v in self.y])
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
        return {nf for f in funcs for nf in f.extract(NamedFunction)}

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
            key=lambda nf: (len(nf.extract(NamedFunction)),
                            nf.name, len(nf.ders), nf.ders)))

    @cached_property
    def c(self) -> tuple[ConstName, ...]:
        return tuple(sorted(
            {c for f in self.a_set for c in f.extract(ConstName)} ))

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
    def _const(self, f: Const) -> str:
        """
        Const as float with 'dp' kind.
        """
        return f'{f.number}_dp'

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
                f'n{len(f.extract(NamedFunction))}',
                '_'.join([f'd{self.y.index(v)}' for v in f.ders]),
                f.name,
            ])

    def _named_function_definition(self, nf: NamedFunction) -> str:
        return _TEMPLATE_FUNCTION_DEFINITION.format(
                name=self._get_function_name(nf),
                body=
                    f'  res = {self.format(nf.f)}' \
                        if isinstance(nf, ExplicitFunction) else
                        '  ! missing opaque definition')

    def _make_update_a(self) -> str:
        assignments = [f'  {self.format(nf)} = {self._get_function_name(nf)}(y)' \
            for nf in self.a_eval_order]
        return _TEMPLATE_UPDATE_A.format('\n'.join(assignments))

    def _make_f(self) -> str:
        assignments = [f'  ydot({1+i}) = {self.format(f)}' for i, f in enumerate(self.f)]
        return _TEMPLATE_F.format('\n'.join(assignments))

    def _make_jac(self) -> str:
        assignments = []
        for i, row in enumerate(self.jacobian):
            for j, f in enumerate(row):
                if f != ZERO:
                    assignments.append(f'  pd({1+i},{1+j}) = {self.format(f)}')
        return _TEMPLATE_JAC.format('\n'.join(assignments))

    def _make_function_definitions(self) -> str:
        return '\n\n'.join([self._named_function_definition(cf) for cf in self.a_human_order])

    def make_functions_file(self) -> str:
        return _TEMPLATE_FUNCTIONS_FILE.format(
            f=self._make_f(),
            jac=self._make_jac(),
            update_a=self._make_update_a(),
            functions=self._make_function_definitions(),
            )

    def make_program_file(self,
                          *,
                          y0: Iterable[float],
                          t0: float = 0.,
                          consts: dict[ConstName, float],
                          tout: float,
                          rtol: float | Iterable[float],
                          atol: float | Iterable[float],
                          num_steps: int = 1,
                          tout_multiplier: float = 10.,
                          mf: int = 21,
                          itask: int = 1,
                          iopt: int = 0,
                          ) -> str:

        return lsode_program.make_program_file(
                  neq=self.neq,
                  n_yca=self.n_yca,
                  y0=y0,
                  t0=t0,
                  consts=((self.neq + i, c.name, consts[c]) for i, c in enumerate(self.c)),
                  tout=tout,
                  rtol=rtol,
                  atol=atol,
                  num_steps=num_steps,
                  tout_multiplier=tout_multiplier,
                  mf=mf,
                  itask=itask,
                  iopt=iopt,
                  )


def _is_trivial(f: Function) -> bool:
    return isinstance(f, Const | ConstName | Var)

def _unpack_if_trivial(f: Function) -> Function:
    if isinstance(f, ExplicitFunction) and _is_trivial(f.f):
        return f.f
    return f

_TEMPLATE_FUNCTION_DEFINITION = """\
real(dp) function {name}(y) result(res)
  real(dp), intent(in) :: y(*)
{body}
end function"""

_TEMPLATE_UPDATE_A = """\
subroutine update_a(y)
  real(dp), intent(inout) :: y(*)
{}
end subroutine"""

_TEMPLATE_F = """\
subroutine f(neq, t, y, ydot)
  integer, intent(in) :: neq
  real(dp), intent(in) :: t
  real(dp), intent(inout) :: y(*)
  real(dp), intent(out) :: ydot(*)
  call update_a(y)
{}
end subroutine"""

_TEMPLATE_JAC = """\
subroutine jac(neq, t, y, ml, mu, pd, nrpd)
  integer, intent(in) :: neq, ml, mu, nrpd
  real(dp), intent(in) :: t, y(*)
  real(dp), intent(out) :: pd(nrpd, *)
{}
end subroutine"""

_TEMPLATE_FUNCTIONS_FILE = """\
module functions
  use, intrinsic :: iso_c_binding, only: dp => c_double
  implicit none
contains

{f}

{jac}

{update_a}

{functions}

end"""
