from __future__ import annotations
from functools import singledispatchmethod

from formatter_base import Formatter
from chain import (
    Function, Const, ConstName, Var, Sum, Prod, Times, Minus, Frac, Pow, Sqrt,
    Inv, Sin, Cos, Arctan, NamedFunction, ExplicitFunction, OpaqueFunction
)

class FortranFormatter(Formatter):
    """Format Functions into simple Fortran-style code/expression.

    This is intentionally conservative and aims for readable, valid
    Fortran-like expressions (real signatures using real(8)).
    """

    @property
    def prefix_der(self) -> str: return self._prefix_der

    def __init__(self, prefix_der: str = "_d") -> None:
        self._prefix_der = prefix_der

    # Repeating needed because @singledispatchmethod doesn't work well with inheritance
    @singledispatchmethod
    def format(self, f: Function) -> str:
        return super().format(f)

    @format.register
    def _const(self, f: Const) -> str:
        return f'{f.number}D0'

    @format.register
    def _const_name(self, f: ConstName) -> str:
        return str(f.name)

    @format.register
    def _var(self, f: Var) -> str:
        return f.name

    @format.register
    def _sum(self, f: Sum) -> str:
        parts = [self.format(x) for x in f.args]
        s = ' + '.join(parts)
        s = s.replace('+ -', '- ')
        return f'({s})'

    @format.register
    def _prod(self, f: Prod) -> str:
        return '*'.join([self.format(x) for x in f.args])

    @format.register
    def _times(self, f: Times) -> str:
        return f"{self.format(f.n)}*{self.format(f.f)}"

    @format.register
    def _minus(self, f: Minus) -> str:
        return f"-{self.format(f.f)}"

    @format.register
    def _frac(self, f: Frac) -> str:
        return f"({self.format(f.num)})/({self.format(f.den)})"

    @format.register
    def _pow(self, f: Pow) -> str:
        # Fortran uses ** for power
        if f.exp.number.is_integer():
            return f"({self.format(f.base)})**{int(f.exp.number)}"
        return f"({self.format(f.base)})**{self.format(f.exp)}"

    @format.register
    def _sqrt(self, f: Sqrt) -> str:
        return f"sqrt({self.format(f.f)})"

    @format.register
    def _inv(self, f: Inv) -> str:
        return f"1/({self.format(f.f)})"

    @format.register
    def _sin(self, f: Sin) -> str:
        return f"sin({self.format(f.f)})"

    @format.register
    def _cos(self, f: Cos) -> str:
        return f"cos({self.format(f.f)})"

    @format.register
    def _atan(self, f: Arctan) -> str:
        return f"atan({self.format(f.f)})"

    def _get_function_str_with_ders(self, f: NamedFunction) -> str:
        return self.prefix_der.join([f.name, *[str(v) for v in f.ders]])

    @format.register
    def _explicit_function(self, f: ExplicitFunction) -> str:
        return self._get_function_str_with_ders(f)

    @format.register
    def _opaque_function(self, f: OpaqueFunction) -> str:
        return self._get_function_str_with_ders(f)
