from __future__ import annotations
from typing import Set

from formatter_base import Formatter
from chain import (
    Function, Const, Var, Sum, Prod, Times, Minus, Frac, Pow,
    Sq, Sqrt, Inv, Sin, Cos, Arctan, ExplicitFunction, OpaqueFunction
)


def _collect_vars(f: Function, seen: Set[str] | None = None) -> Set[str]:
    if seen is None:
        seen = set()
    if isinstance(f, Var):
        seen.add(f.name)
    # recursively traverse attributes we know about
    if hasattr(f, 'args'):
        for a in getattr(f, 'args'):
            _collect_vars(a, seen)
    if hasattr(f, 'f'):
        _collect_vars(getattr(f, 'f'), seen)
    if hasattr(f, 'num'):
        _collect_vars(getattr(f, 'num'), seen)
    if hasattr(f, 'den'):
        _collect_vars(getattr(f, 'den'), seen)
    return seen


class FortranFormatter(Formatter):
    """Format Functions into simple Fortran-style code/expression.

    This is intentionally conservative and aims for readable, valid
    Fortran-like expressions (real signatures using real(8)).
    """

    # Primitive types
    @Formatter.format.register
    def _const(self, f: Const) -> str:
        return str(f.number)

    @Formatter.format.register
    def _var(self, f: Var) -> str:
        return f.name

    # Composite
    @Formatter.format.register
    def _sum(self, f: Sum) -> str:
        parts = [self.format(x) for x in f.args]
        s = ' + '.join(parts)
        s = s.replace('+ -', '- ')
        return s

    @Formatter.format.register
    def _prod(self, f: Prod) -> str:
        return '*'.join([self.format(x) for x in f.args])

    @Formatter.format.register
    def _times(self, f: Times) -> str:
        # n*expr or constant folding
        expr = self.format(f.f)
        return f"{self.format(f.n)}*{expr}"

    @Formatter.format.register
    def _minus(self, f: Minus) -> str:
        return f"-{self.format(f.f)}"

    @Formatter.format.register
    def _frac(self, f: Frac) -> str:
        return f"({self.format(f.num)})/({self.format(f.den)})"

    @Formatter.format.register
    def _pow(self, f: Pow) -> str:
        # Fortran uses ** for power
        return f"({f.base})**{f.exp}"

    @Formatter.format.register
    def _sqrt(self, f: Sqrt) -> str:
        return f"sqrt({self.format(f.f)})"

    @Formatter.format.register
    def _inv(self, f: Inv) -> str:
        return f"1/({self.format(f.f)})"

    @Formatter.format.register
    def _sin(self, f: Sin) -> str:
        return f"sin({self.format(f.f)})"

    @Formatter.format.register
    def _cos(self, f: Cos) -> str:
        return f"cos({self.format(f.f)})"

    @Formatter.format.register
    def _atan(self, f: Arctan) -> str:
        return f"atan({self.format(f.f)})"

    # Explicit / Opaque functions: format as call, but also provide definitions
    @Formatter.format.register
    def _explicit_function(self, f: ExplicitFunction) -> str:
        args = _collect_vars(f.f)
        args_list = ', '.join(sorted(args))
        return f"{f.name}({args_list})"

    @Formatter.format.register
    def _opaque_function(self, f: OpaqueFunction) -> str:
        args_list = ', '.join([v.name for v in f.vars_])
        return f"{f.name}({args_list})"

    def explicit_function_definition(self, ef: ExplicitFunction) -> str:
        vars_ = sorted(list(_collect_vars(ef.f)))
        args_sig = ', '.join(vars_)
        lines = []
        lines.append(f"function {ef.name}({args_sig}) result(res)")
        lines.append("  real(8) :: res")
        if vars_:
            decl = ', '.join(vars_)
            lines.append(f"  real(8), intent(in) :: {decl}")
        lines.append(f"  res = {self.format(ef.f)}")
        lines.append(f"end function {ef.name}")
        return '\n'.join(lines)

    def opaque_function_definition(self, of: OpaqueFunction) -> str:
        args = ', '.join([v.name for v in of.vars_])
        lines = [f"! Opaque function {of.name}({args}) -- signature only"]
        lines.append(f"function {of.name}({args}) result(res)")
        lines.append("  real(8) :: res")
        if args:
            lines.append(f"  real(8), intent(in) :: {args}")
        lines.append(f"  ! Implementation not provided (opaque)")
        lines.append(f"end function {of.name}")
        return '\n'.join(lines)
