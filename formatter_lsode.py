from __future__ import annotations
from functools import singledispatchmethod, cached_property

from formatter_fortran import FortranFormatter
from chain import (
    Function, ConstName, Var, NamedFunction, ExplicitFunction, OpaqueFunction,
    extract
)

class LsodeFormatter(FortranFormatter):
    """
    Takes FortranFormatter further to allow LSODE use.

    Introduces the concept of variable vector (y), function matrix (F),
    Jacobian matrix and implements function definitions

    solution_vector: the "y" of LSODE
    function_vector: the "f" of y'=f(y)
    constant_vector: extracted ConstNames go here for reference in Fortran
    named_functions (aux vector): extracted NamedFunctions go here for memoizing and reference in Fortran
    """

    @property
    def solution_vector(self) -> list[Var]: return self._solution_vector

    @cached_property
    def function_vector(self) -> list[Function]:
        l = []
        for i, f in enumerate(self._function_vector):
            l.append(f if isinstance(f, Var) else
                        ExplicitFunction(f'{self.function_vector_name}{1+i}', f))
        return l

    @property
    def prefix_function(self) -> str: return self._prefix_function

    @property
    def solution_vector_name(self) -> str: return self._solution_vector_name

    @property
    def function_vector_name(self) -> str: return self._function_vector_name

    @property
    def constant_vector_name(self) -> str: return self._constant_vector_name

    @property
    def aux_vector_name(self) -> str: return self._aux_vector_name

    @property
    def float_type(self) -> str: return self._float_type

    def __init__(self,
            solution_vector: list[Var],
            function_vector: list[Var],
            *,
            prefix_der: str = "_d",
            prefix_function: str = "__",
            solution_vector_name: str = "y",
            function_vector_name: str = "f",
            constant_vector_name: str = "c",
            aux_vector_name: str = "a",
            float_type: str = "double precision",
            ) -> None:
        super().__init__(prefix_der)
        self._solution_vector = solution_vector
        self._function_vector = function_vector
        self._prefix_function = prefix_function
        self._solution_vector_name = solution_vector_name
        self._function_vector_name = function_vector_name
        self._constant_vector_name = constant_vector_name
        self._aux_vector_name = aux_vector_name
        self._float_type = float_type

    # Shorthands

    @property
    def y_vec(self) -> list[Var]: return self.solution_vector

    @property
    def y_str(self) -> str: return self.solution_vector_name

    @property
    def f_vec(self) -> list[Var]: return self.function_vector

    @property
    def f_str(self) -> str: return self.function_vector_name

    @property
    def c_vec(self) -> list[Var]: return self.constant_vector

    @property
    def c_str(self) -> str: return self.constant_vector_name

    @property
    def aux_vec(self) -> list[Var]: return self.named_functions_eval_order

    @property
    def aux_str(self) -> str: return self.aux_vector_name

    @cached_property
    def jacobian(self) -> list[list[Function]]:
        return [[f.der(v) for v in self.y_vec] for f in self.f_vec]

    @cached_property
    def named_functions(self) -> set[NamedFunction]:
        """
        Visits everyone in function vector + jacobian to find all NamedFunctions used.

        The result is sorted so that opaque functions come first, for ease of
        finding them to be manually defined later.
        """
        funcs = set(self.f_vec)
        funcs.update({f for row in self.jacobian for f in row})
        return {nf for f in funcs for nf in extract(f, NamedFunction)}

    @property
    def named_functions_human_order(self) -> tuple[NamedFunction, ...]:
        """
        NamedFunction objects in alphabetical order (but opaque and function
        vector first), for ease of manually defining.
        """
        return tuple(sorted(self.named_functions,
            key=lambda nf: (1 if isinstance(nf, OpaqueFunction) \
                            else 2 if nf in self.f_vec \
                            else 3,
                            nf.name, len(nf.ders), nf.ders)))

    @cached_property
    def named_functions_eval_order(self) -> tuple[NamedFunction, ...]:
        """
        NamedFunction objects in the correct order of evaluation.

        NamedFunction in the optimized order of evaluation, which allows
        saving each function in a vector to be used by the next functions. That
        is, each is the building block of the next. It is done by counting how
        many functions each uses, and ordering ascending.
        """
        return tuple(sorted(self.named_functions,
            key=lambda nf: (len(extract(nf, NamedFunction)),
                            nf.name, len(nf.ders), nf.ders)))

    @cached_property
    def constant_vector(self) -> list[ConstName]:
        return sorted({c for f in self.named_functions for c in extract(f, ConstName)})

    # Repeating needed because @singledispatchmethod doesn't work well with inheritance
    @singledispatchmethod
    def format(self, f: Function) -> str:
        return super().format(f)

    # Overrides for using vectors, e.g turning Var "XG" into "y(1)"

    @format.register
    def _var(self, f: Var) -> str:
        return f'{self.y_str}({1+self.y_vec.index(f)})'

    @format.register
    def _const_name(self, f: ConstName) -> str:
        return f'{self.c_str}({1+self.c_vec.index(f)})'

    @format.register
    def _named(self, f: NamedFunction) -> str:
        return f'{self.aux_str}({1+self.aux_vec.index(f)})'

    # Repeating these two is needed because else the more specific ones from
    # FortranFormatter are called

    @format.register
    def _named_explicit(self, f: ExplicitFunction) -> str:
        return self._named(f)

    @format.register
    def _named_opaque(self, f: OpaqueFunction) -> str:
        return self._named(f)

    def _get_function_name(self, f: NamedFunction) -> str:
        return self.prefix_function + self._get_function_str_with_ders(f)

    def _get_function_args(self):
        return f'({self.y_str}, {self.c_str}, {self.aux_str})'

    def get_function_call(self, f: NamedFunction) -> str:
        return self._get_function_name(f) + self._get_function_args()

    def _get_function_header(self, f: NamedFunction) -> str:
        lines = []
        lines.append(f'function {self.get_function_call(f)} result(res)')
        lines.append(f'  {self.float_type} :: res')
        lines.append(f'  {self.float_type}, intent(in) :: {self.y_str}({len(self.y_vec)})')
        lines.append(f'  {self.float_type}, intent(in) :: {self.c_str}({len(self.c_vec)})')
        lines.append(f'  {self.float_type}, intent(in) :: {self.aux_str}({len(self.aux_vec)})')
        return '\n'.join(lines)

    def _get_entire_function(self, f: NamedFunction, body: str) -> str:
        lines = []
        lines.append(body)
        lines.append('end function')
        return '\n\n'.join([self._get_function_header(f), '\n'.join(lines)])

    def named_function_definition(self, nf: NamedFunction) -> str:
        body = ('  res = ' + self.format(nf.f)) \
            if isinstance(nf, ExplicitFunction) else '  ! missing opaque definition'
        return self._get_entire_function(nf, body)

    def make_update_aux(self):
        lines = []
        lines.append(f'subroutine update_aux{self._get_function_args()}')
        lines.append(f'  {self.float_type}, intent(in) :: {self.y_str}({len(self.y_vec)})')
        lines.append(f'  {self.float_type}, intent(in) :: {self.c_str}({len(self.c_vec)})')
        lines.append(f'  {self.float_type}, intent(inout) :: {self.aux_str}({len(self.aux_vec)})')
        lines.append('')
        for i, nf in enumerate(self.named_functions_eval_order):
            lines.append(f'  {self.aux_str}[{1+i}] = {self.get_function_call(nf)}')
        lines.append('end subroutine')
        return '\n'.join(lines)
