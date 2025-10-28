from __future__ import annotations
from functools import cached_property

from formatter_fortran import FortranFormatter
from chain import (
    Function, Var, NamedFunction, ExplicitFunction, OpaqueFunction,
    extract
)

class LsodeFormatter(FortranFormatter):
    """Takes FortranFormatter further to allow LSODE use.

    Introduces the concept of variable vector (y), function matrix (F),
    Jacobian matrix and implements function definitions
    """

    @property
    def solution_vector(self) -> list[Var]: return self._solution_vector

    @property
    def function_vector(self) -> list[Var]: return self._function_vector

    @property
    def prefix_function(self) -> str: return self._prefix_function

    @property
    def solution_vector_name(self) -> str: return self._solution_vector_name

    @property
    def function_vector_name(self) -> str: return self._function_vector_name

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
            float_type: str = "double precision",
            ) -> None:
        super().__init__(prefix_der)
        self._solution_vector = solution_vector
        self._function_vector = function_vector
        self._prefix_function = prefix_function
        self._solution_vector_name = solution_vector_name
        self._function_vector_name = function_vector_name
        self._float_type = float_type

    @property
    def y_vec(self) -> list[Var]: return self.solution_vector

    @property
    def y_str(self) -> str: return self.solution_vector_name

    @property
    def n(self) -> int: return len(self.solution_vector)

    @property
    def f_vec(self) -> list[Var]: return self.function_vector

    @property
    def f_str(self) -> str: return self.function_vector_name

    @cached_property
    def jacobian(self) -> list[list[Function]]:
        return [[f.der(v) for v in self.y_vec] for f in self.f_vec]

    @cached_property
    def named_functions(self) -> list[NamedFunction]:
        funcs = set(self.f_vec)
        funcs.update({f for row in self.jacobian for f in row})
        return sorted(
                {nf for f in funcs for nf in extract(f, NamedFunction)},
                key=_named_functions_key)

    @FortranFormatter.format.register
    def _var(self, f: Var) -> str:
        """
        Custom implementation for Var: Converts Var('x') to 'Y(1)'.
        """
        return f'{self.y_str}({1+self.y_vec.index(f)})'

    def _get_function_name(self, f: NamedFunction) -> str:
        return self.prefix_function + self._get_function_str_with_ders(f)

    def _get_function_call(self, f: NamedFunction) -> str:
        return self._get_function_name(f) + f'({self.y_str})'

    @FortranFormatter.format.register
    def _explicit_function(self, f: ExplicitFunction) -> str:
        return self._get_function_call(f)

    @FortranFormatter.format.register
    def _opaque_function(self, f: OpaqueFunction) -> str:
        return self._get_function_call(f)

    def _get_function_header(self, f: NamedFunction) -> str:
        lines = []
        lines.append(f'function {self._get_function_name(f)}({self.y_str}) result(res)')
        lines.append(f'  {self.float_type} :: res')
        lines.append(f'  {self.float_type}, intent(in) :: {self.y_str}({self.n})')
        return '\n'.join(lines)

    def _get_entire_function(self, f: NamedFunction, body: str) -> str:
        lines = []
        lines.append(body)
        lines.append('end function')
        return '\n\n'.join([self._get_function_header(f), '\n'.join(lines)])

    def named_function_definition(self, nf: NamedFunction) -> str:
        if isinstance(nf, ExplicitFunction):
            return self._get_entire_function(nf, '  res = ' + self.format(nf.f))
        elif isinstance(nf, OpaqueFunction):
            return self._get_entire_function(nf, '  ! missing opaque definition')
        else:
            raise TypeError()

def _named_functions_key(nf: NamedFunction) -> list:
    if isinstance(nf, OpaqueFunction):
        return (1, nf.name, len(nf.ders), nf.ders)
    elif isinstance(nf, ExplicitFunction):
        return (2, nf.name, len(nf.ders), nf.ders)
    else:
        raise TypeError()
