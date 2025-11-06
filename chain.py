from __future__ import annotations
# this import is so Function can reference itself

from abc import ABC, abstractmethod
from typing import cast, Any, Iterable, Callable
from collections.abc import Hashable
from functools import cache
from dataclasses import dataclass
import math
import dataclasses

class _HasName:
    @property
    def name(self) -> str: return self._name
    def __init__(self, name: str) -> None: self._name = name

class _HasF:
    @property
    def f(self) -> Function: return self._f
    def __init__(self, f: Function) -> None: self._f = f

class _HasN:
    @property
    def n(self) -> Const: return self._n
    def __init__(self, n: Const) -> None: self._n = n

class _HasVars:
    @property
    def vars_(self) -> tuple[Var, ...]: return self._vars_
    def __init__(self, vars_: tuple[Var, ...]) -> None: self._vars_ = vars_

class _HasDers:
    @property
    def ders(self) -> tuple[Var, ...]: return self._ders
    def __init__(self, ders: tuple[Var, ...]) -> None: self._ders = ders

class _HasKey(ABC):
    @abstractmethod
    def key(self) -> Any:
        pass

class _ExtractableMixin(_HasKey):
    # Hashable is needed for the @cache decorator in _extract_any, but this is a
    # bug (mypy bugs 11469 and 11470), since type objects are always hashable
    def extract[T: Hashable](self, t: type[T]) -> set[T]:
        """
        Recursively extracts all instances of a given type found in self and
        children (self.key).  This is key for formatters (like LsodeFormatter)
        to do important things like find all ConstNames of a system or build
        function dependency graphs.
        """
        s1 = {self} if isinstance(self, t) else set()
        s2 = _extract_any(self.key(), t)
        return s1 | s2

class Function(_ExtractableMixin, _HasKey, ABC):
    @abstractmethod
    def der(self, v: Var) -> Function:
        pass

    def __lt__(self, other: object) -> bool:
        if isinstance(other, float | int):
            other = _k(other)
        elif not isinstance(other, Function):
            return NotImplemented
        t1 = type(self)
        t2 = type(other)
        if t1 != t2:
            return _SORT_ORDER.index(t1) < _SORT_ORDER.index(t2)
        else:
            return self.key() < other.key() # type: ignore

    def __eq__(self, other: object) -> bool:
        if isinstance(other, float | int):
            other = _k(other)
        if not isinstance(other, Function):
            return False
        return type(self) == type(other) and self.key() == other.key()

    def __gt__(self, other: object) -> bool:
        return not self.__lt__(other) and not self.__eq__(other)
    def __le__(self, other: object) -> bool:
        return self.__lt__(other) or self.__eq__(other)
    def __ge__(self, other: object) -> bool:
        return not self.__lt__(other)
    def __ne__(self, other: object) -> bool:
        return not self.__eq__(other)

    def __add__(self, other: object) -> Function:
        if not isinstance(other, Function | float | int): return NotImplemented
        return Sum.factory(self, other)

    def __mul__(self, other: object) -> Function:
        if not isinstance(other, Function | float | int): return NotImplemented
        return Prod.factory(self, other)

    def __truediv__(self, other: object) -> Function:
        if not isinstance(other, Function | float | int): return NotImplemented
        return Frac.factory(self, other)

    def __rtruediv__(self, other: object) -> Function:
        if not isinstance(other, Function | float | int): return NotImplemented
        return Frac.factory(other, self)

    def __pow__(self, other: object) -> Function:
        if not isinstance(other, Const | float | int): return NotImplemented
        return Pow.factory(self, other)

    def __sub__(self, other: object) -> Function:
        if not isinstance(other, Function | float | int): return NotImplemented
        return Sum.factory(self, Minus.factory(other))

    def __rsub__(self, other: object) -> Function:
        if not isinstance(other, Function | float | int): return NotImplemented
        return Sum.factory(other, Minus.factory(self))

    def __neg__(self) -> Function:
        return Minus.factory(self)

    __radd__ = __add__
    __rmul__ = __mul__

    def __repr__(self) -> str:
        key = repr(self.key())
        if key[0] != '(':
            key = f'({key})'
        return f'{type(self).__name__}{key}'

    def __hash__(self) -> int:
        return hash(( type(self), self.key() ))

class FunctionVariadic(Function):
    @property
    def args(self) -> tuple[Function, ...]:
        return self._args

    def key(self) -> tuple[Function, ...]:
        return self.args

    def __init__(self, *args: Function) -> None:
        self._args = args

    @staticmethod
    @abstractmethod
    def factory(*args: Function) -> Function:
        pass

class FunctionWithSingleFunction(Function, _HasF):
    def __init__(self, f: Function) -> None:
        _HasF.__init__(self, f)

    @abstractmethod
    def copy(self, f: Function) -> Function:
        """Reconstruct with another function argument"""
        pass

    @property
    @abstractmethod
    def selfder(self) -> Function:
        """Chain Rule"""
        pass

    def der(self, v: Var) -> Function:
        return Prod.factory(self.selfder, self.f.der(v))

class FunctionSimple(FunctionWithSingleFunction):
    @staticmethod
    @abstractmethod
    def factory(f: Function) -> Function:
        pass

    def copy(self, f: Function) -> Function:
        return self.factory(f)

    def key(self) -> Function:
        return self.f

class FunctionWithCardinality(FunctionWithSingleFunction, _HasN):
    def key(self) -> tuple[Function, ...]:
        return (self.f, self.n)

    def __init__(self, f: Function, n: Const):
        _HasF.__init__(self, f)
        _HasN.__init__(self, n)

    def copy(self, f: Function) -> Function:
        return self.factory(f, self.n)

    @staticmethod
    @abstractmethod
    def factory(f: Function, n: Const) -> Function:
        pass

class Const(Function):
    @property
    def number(self) -> float:
        return self._number

    @staticmethod
    def factory(number: float) -> Const:
        if number is True:
            number = 1
        elif number is False:
            number = 0

        return Const(number)

    def __init__(self, number: float) -> None:
        self._number = float(number)

    def key(self) -> float:
        return self.number

    def der(self, v: Var) -> Const:
        return ZERO

    def __str__(self) -> str:
        if self.number.is_integer():
            return str(int(self.number))
        else:
            return str(self.number)

    def __add__(self, other: object) -> Const:
        if not isinstance(other, Const | float | int): return NotImplemented
        return _k(self.number + _fix_const(other).number)

    def __mul__(self, other: object) -> Const:
        if not isinstance(other, Const | float | int): return NotImplemented
        return _k(self.number * _fix_const(other).number)

    def __truediv__(self, other: object) -> Const:
        if not isinstance(other, Const | float | int): return NotImplemented
        return _k(self.number / _fix_const(other).number)

    def __rtruediv__(self, other: object) -> Const:
        if not isinstance(other, Const | float | int): return NotImplemented
        return _k(_fix_const(other).number / self.number)

    def __pow__(self, other: object) -> Const:
        if not isinstance(other, Const | float | int): return NotImplemented
        return _k(self.number ** _fix_const(other).number)

    def __rpow__(self, other: object) -> Const:
        if not isinstance(other, Const | float | int): return NotImplemented
        return _k(_fix_const(other).number ** self.number)

    def __sub__(self, other: object) -> Const:
        if not isinstance(other, Const | float | int): return NotImplemented
        return _k(self.number - _fix_const(other).number)

    def __rsub__(self, other: object) -> Const:
        if not isinstance(other, Const | float | int): return NotImplemented
        return _k(_fix_const(other).number - self.number)

    def __neg__(self) -> Const:
        return _k(-self.number)

    __radd__ = __add__
    __rmul__ = __mul__

ZERO     : Const = Const.factory(0)
ONE      : Const = Const.factory(1)
TWO      : Const = Const.factory(2)
HALF     : Const = Const.factory(.5)
MINUS_ONE: Const = Const.factory(-1)
PI       : Const = Const.factory(math.pi)

class ConstName(Function, _HasName):
    def __init__(self, name: str) -> None:
        _HasName.__init__(self, name)

    def key(self) -> str:
        return self.name

    def der(self, v: Var) -> Const:
        return ZERO

    def __str__(self) -> str:
        return self.name

class Var(Function, _HasName):
    def __init__(self, name: str) -> None:
        _HasName.__init__(self, name)

    def key(self) -> str:
        return self.name

    def der(self, v: Var) -> Const:
        return _k(v == self)

    def __str__(self) -> str:
        return self.name

X: Var = Var('x')
Y: Var = Var('y')
Z: Var = Var('z')

class Sum(FunctionVariadic):
    @staticmethod
    def factory(*args: Function | float) -> Function:
        if not args:
            return ZERO

        args_list = _fix_consts_in_iterable(args)
        del args

        if all([isinstance(x, Const) for x in args_list]):
            return _k(math.fsum([cast(Const, c).number for c in args_list]))

        lc: list[Const]
        lc, lf = dataclasses.astuple(_flatten_sum(*args_list))
        k = math.fsum([c.number for c in lc])
        args_list = lf if k == 0 else [_k(k), *lf]
        del k, lf

        args_list = _fold(args_list, Sum, Times)

        if len(args_list) == 0:
            return ZERO
        elif len(args_list) == 1:
            return args_list[0]
        else:
            return Sum(*sorted(args_list))

    def der(self, v: Var) -> Function:
        return Sum.factory(*[x.der(v) for x in self.args])

    def __str__(self) -> str:
        return ' + '.join([str(x) for x in self.args]).replace(' + -', ' - ')

class Prod(FunctionVariadic):
    @staticmethod
    def factory(*args: Function | float) -> Function:
        if not args:
            return ONE

        args_list = _fix_consts_in_iterable(args)
        del args

        if all([isinstance(x, Const) for x in args_list]):
            return _k(math.prod([cast(Const, c).number for c in args_list]))

        lc: list[Const]
        nums: list[Function]
        dens: list[Function]
        lc, nums, dens = dataclasses.astuple(_flatten_prod(*args_list))
        k = math.prod([c.number for c in lc])
        if k == 0:
            return ZERO
        del args_list, lc

        if dens:
            den = Prod.factory(*dens)
            if den != ONE:
                num = Prod.factory(*nums)
                return Times.factory(Frac.factory(num, den), k)
            del den
        del dens

        nums = _fold(nums, Prod, Pow)

        if len(nums) == 0:
            return _k(k)
        elif len(nums) == 1:
            return Times.factory(nums[0], k)
        else:
            return Times.factory(Prod(*sorted(nums)), k)

    def der(self, v: Var) -> Function:
        list_to_sum: list[Function] = []
        for i, _ in enumerate(self.args):
            list_to_multiply: list[Function] = []
            for j, f in enumerate(self.args):
                if i == j:
                    list_to_multiply.append(f.der(v))
                else:
                    list_to_multiply.append(f)
            list_to_sum.append(Prod.factory(*list_to_multiply))

        return Sum.factory(*list_to_sum)

    def __str__(self) -> str:
        return '*'.join([_parens_if_sum(x) for x in self.args])

class Times(FunctionWithCardinality):
    @staticmethod
    def factory(f: Function, nf: Const | float) -> Function:
        n = _fix_const(nf)
        del nf

        if isinstance(f, Times):
            n = n*f.n
            f = f.f

        if n == ZERO:
            return ZERO
        elif isinstance(f, Const):
            return _k(n.number * f.number)
        elif n == ONE:
            return f
        elif n == MINUS_ONE:
            return Minus.factory(f)

        return Times(f, n)

    @property
    def selfder(self) -> Const:
        return self.n

    def __str__(self) -> str:
        if isinstance(self.f, Sum):
            return f'{self.n}*({self.f})'
        else:
            return f'{self.n}*{self.f}'

class Minus(Times):
    @staticmethod
    def factory(
            ff: Function | float,
            nf: Const | float = MINUS_ONE
            ) -> Function:
        f = _fix_const_if_float(ff)
        n = _fix_const(nf)
        del ff, nf

        if n != MINUS_ONE:
            raise ValueError('Minus cardinality must be -1')

        if isinstance(f, Minus):
            return f.f
        elif isinstance(f, Const):
            return _k(-f.number)
        elif isinstance(f, Times):
            return Times.factory(f.f, f.n.number * -1)

        return Minus(f)

    def __init__(self, f: Function):
        super().__init__(f, MINUS_ONE)

    def __str__(self) -> str:
        if isinstance(self.f, Sum):
            return f'-({self.f})'
        else:
            return f'-{self.f}'

class Frac(Function):
    @staticmethod
    def factory(numf: Function | float, denf: Function | float) -> Function:
        num = _fix_const_if_float(numf)
        den = _fix_const_if_float(denf)
        del numf, denf

        if den == ZERO:
            raise ZeroDivisionError('denominator cannot be zero')
        elif den == ONE:
            return num
        elif isinstance(num, Const):
            if isinstance(den, Const):
                return _k(num.number/den.number)
            else:
                return Times.factory(Inv.factory(den), num)
        elif isinstance(den, Const):
            return Times.factory(num, _k(1/den.number))

        lc: list[Const]
        nums: list[Function]
        dens: list[Function]
        lc, nums, dens = dataclasses.astuple(_flatten_prod(Frac(num, den)))
        k = math.prod([c.number for c in lc])
        if k == 0:
            return ZERO
        del lc

        dca = _fold_make_dca(nums, dens, Pow)
        if not dca:
            return _k(k)
        del nums, dens

        dca_nums: dict[float, Function] = {}
        dca_dens: dict[float, Function] = {}
        for (n, lf) in dca.items():
            f = Prod.factory(*lf)
            if   n > 0: dca_nums[ n] = f
            elif n < 0: dca_dens[-n] = f
        del n, lf, f

        # Warning: from here on, we use direct constructors instead
        # of factories where preventing infinite recursions is needed
        # e.g: Frac(num,den) instead of Frac.factory(num,den)

        fracs: list[Function] = []
        nums = []
        dens = []
        for n, num in dca_nums.items():
            if n in dca_dens:
                den = dca_dens[n]
                fracs.append(Pow.factory(Frac(num, den), n))
            else:
                nums.append(Pow.factory(num, n))
        for n, den in dca_dens.items():
            if n not in dca_nums:
                dens.append(Pow.factory(den, n))
        del dca_nums, dca_dens, num, den

        prod_list: list[Function]
        if not dens:
            prod_list = [*fracs, *nums]
        elif not nums:
            prod_list = [*fracs, Inv.factory(Prod.factory(*dens))]
        else:
            prod_list = [*fracs, Frac(Prod.factory(*nums), Prod.factory(*dens))]
        prod_list = sorted(prod_list)

        if len(prod_list) > 1:
            return Times.factory(Prod(*prod_list), k)
        else:
            return Times.factory(prod_list[0], k)

    def key(self) -> tuple[Function, ...]:
        return (self.num, self.den)

    def __init__(self, num: Function, den: Function):
        self.num = num
        self.den = den

    def der(self, v: Var) -> Function:
        return Frac.factory(
                Sum.factory(
                    Prod.factory(self.num.der(v), self.den),
                    Minus.factory(
                        Prod.factory(self.num, self.den.der(v)),
                )),
                Sq.factory(self.den))

    def __str__(self) -> str:
        den = str(self.den) if isinstance(self.den, Minus) else \
                _parens_if_sum_prod_frac_times(self.den)
        return f'{_parens_if_sum(self.num)}/{den}'

class Pow(FunctionWithCardinality):
    @staticmethod
    def factory(base: Function, expf: Const | float) -> Function:
        exp = _fix_const(expf)
        del expf

        if isinstance(base, Pow):
            exp = _k(base.exp.number * exp.number)
            base = base.base

        if base == ZERO and exp == ZERO:
            raise ValueError('0^0 is undefined')
        elif base == ZERO and exp < ZERO:
            raise ValueError('0 raised to a negative exponent is undefined')
        elif base == ZERO:
            return ZERO
        elif exp == ZERO:
            return ONE
        elif exp == ONE:
            return base
        elif isinstance(base, Const):
            if base.number < 0 and not exp.number.is_integer():
                raise ValueError('Cannot have complex values')
            return _k(pow(base.number, exp.number))
        elif isinstance(base, Times):
            return Times.factory(
                    Pow.factory(base.f, exp),
                    base.n.number ** exp.number
            )
        elif isinstance(base, Frac) and exp < ZERO:
            return Pow.factory(Frac.factory(base.den, base.num), -exp)
        elif exp == MINUS_ONE:
            return Inv(base)
        elif exp == TWO:
            return Sq(base)
        elif exp == HALF:
            return Sqrt(base)

        return Pow(base, exp)

    @property
    def base(self) -> Function:
        return self.f

    @property
    def exp(self) -> Const:
        return self.n

    @property
    def selfder(self) -> Function:
        return Times.factory(
                Pow.factory(self.base, self.exp.number - 1),
                self.exp)

    def der(self, v: Var) -> Function:
        # small optimization for sin^2 and cos^2 (derivative = sin(2x))
        base = self.base
        if self.exp == TWO:
            is_sin, is_cos = isinstance(base, Sin), isinstance(base, Cos)
            if is_sin or is_cos:
                # mypy doesn't accept the check above with variables
                f = cast(FunctionSimple, base).f
                ret = Prod.factory(
                        Sin.factory(Times.factory(f, TWO)),
                        f.der(v))
                return ret if is_sin else Minus.factory(ret)

        return super().der(v)

    def __str__(self) -> str:
        return f'{_parens_if_sum_prod_frac_times(self.base)}^{self.exp}'

class Sq(Pow):
    @staticmethod
    def factory(base: Function, expf: Const | float = TWO) -> Function:
        exp = _fix_const(expf)
        del expf

        if exp != TWO:
            raise ValueError('Sq cardinality must be 2')

        return Pow.factory(base, exp)

    def __init__(self, base: Function) -> None:
        super().__init__(base, TWO)

class Sqrt(Pow):
    @staticmethod
    def factory(base: Function, expf: Const | float = HALF) -> Function:
        exp = _fix_const(expf)
        del expf

        if exp != HALF:
            raise ValueError('Sqrt cardinality must be .5')

        return Pow.factory(base, exp)

    def __init__(self, base: Function) -> None:
        super().__init__(base, HALF)

class Inv(Pow):
    @staticmethod
    def factory(f: Function, expf: Const | float = MINUS_ONE) -> Function:
        exp = _fix_const(expf)
        del expf

        if exp != MINUS_ONE:
            raise ValueError('Inv cardinality must be -1')

        return Pow.factory(f, exp)

    def __init__(self, f: Function) -> None:
        super().__init__(f, MINUS_ONE)

    def __str__(self) -> str:
        return f'1/{_parens_if_sum_prod_frac_times(self.f)}'

class Sin(FunctionSimple):
    @staticmethod
    def factory(f: Function) -> Function:
        if isinstance(f, Const):
            if f == ZERO or f == PI or (f.number/math.pi).is_integer():
                return ZERO
            x = 2*f.number/math.pi
            if x.is_integer():
                return ONE if x % 4 == 1 else MINUS_ONE
            return _k(math.sin(f.number))
        return Sin(f)

    @property
    def selfder(self) -> Function:
        return Cos.factory(self.f)

    def __str__(self) -> str:
        return f'sin({self.f})'

class Cos(FunctionSimple):
    @staticmethod
    def factory(f: Function) -> Function:
        if isinstance(f, Const):
            if f == ZERO:
                return ONE
            elif f == PI:
                return MINUS_ONE
            x = 2*f.number/math.pi
            if x.is_integer():
                match x % 4:
                    case 0:
                        return ONE
                    case 2:
                        return MINUS_ONE
                    case 1 | 3:
                        return ZERO
            return _k(math.cos(f.number))
        return Cos(f)

    @property
    def selfder(self) -> Function:
        return Minus.factory(Sin.factory(self.f))

    def __str__(self) -> str:
        return f'cos({self.f})'

class Arctan(FunctionSimple):
    @staticmethod
    def factory(f: Function) -> Function:
        if isinstance(f, Const):
            if f == ZERO:
                return ZERO
            elif f == ONE:
                return _k(math.pi/4)
            elif f == MINUS_ONE:
                return _k(-math.pi/4)
            return _k(math.atan(f.number))
        return Arctan(f)

    @property
    def selfder(self) -> Function:
        return Inv.factory(Sum.factory(ONE, Sq.factory(self.f)))

    def __str__(self) -> str:
        return f'atan({self.f})'

class NamedFunction(Function, _HasName, _HasDers):
    def __init__(self, name: str, ders: tuple[Var, ...] = ()) -> None:
        _HasName.__init__(self, name)
        _HasDers.__init__(self, ders)

class ExplicitFunction(NamedFunction, _HasF):
    def __init__(self, name: str, f: Function, ders: tuple[Var, ...] = ()) -> None:
        super().__init__(name, ders)
        _HasF.__init__(self, f)

    def der(self, v: Var) -> Function:
        d = self.f.der(v)
        return d if d == ZERO else \
            ExplicitFunction(self.name, d, _append_var(self.ders, v))

    def key(self) -> tuple[Any, ...]:
        return (self.name, self.f, self.ders)

    def __str__(self) -> str:
        return _user_function_str(self.name, self.ders)

class OpaqueFunction(NamedFunction, _HasVars):
    def __init__(self, name: str,
            vars_: tuple[Var, ...],
            ders : tuple[Var, ...] = ()) -> None:
        super().__init__(name, ders)
        _HasVars.__init__(self, vars_)

    def der(self, v: Var) -> Function:
        if v in self.vars_:
            return OpaqueFunction(
                self.name, self.vars_, _append_var(self.ders, v))
        else:
            return ZERO

    def key(self) -> tuple[Any, ...]:
        return (self.name, self.vars_, self.ders)

    def __str__(self) -> str:
        return _user_function_str(self.name, self.ders)

def _user_function_str(name: str, ders: tuple[Var, ...]) -> str:
    if not ders:
        #return f'f[{name}]'
        return name
    sders = ', '.join([str(v) for v in ders])
    return f'd[{name}, {sders}]'

def _append_var(tup: tuple[Var, ...], v: Var) -> tuple[Var, ...]:
    return tuple([*tup, v])

def _reduce_consts(input_list: Iterable[Function],
                   func_reduce: Callable[[Iterable[float]], float],
        ) -> tuple[float, list[Function]]:
    consts: list[Const] = []
    rest: list[Function] = []
    for x in input_list:
        if isinstance(x, Const):
            consts.append(x)
        else:
            rest.append(x)
    return (func_reduce([k.number for k in consts]), rest)

# We flatten recursive Sums (or Prods, that is, the passed FunctionVariadic),
# but we also look for its cardinal counterpart (Times for Sum and Pow for Prod)
# to unpack it as well in case it contains a Sum (or Prod).
def _flatten_old(input_list: Iterable[Function],
             variadic_type: type[FunctionVariadic],
             cardinal_type: type[FunctionWithCardinality],
            ) -> list[Function]:
    ret = list(input_list)
    del input_list
    while True:
        must_repeat = False
        new: list[Function] = []
        for x in ret:
            if isinstance(x, variadic_type):
                must_repeat = True
                new.extend(x.args)
            elif isinstance(x, cardinal_type) and isinstance(x.f, variadic_type):
                must_repeat = True
                for y in x.f.args:
                    new.append(cardinal_type.factory(y, x.n))
            else:
                new.append(x)
        ret = new
        if not must_repeat: break
    return ret

# Makes an inverted "args cardinality dict" (ACD -> DCA) from
# the number of appearances of functions in an input list.
# If the passed "cardinal type" appears in the list (e.g: Times
# for a list that is being summed), its n attribute is used.
# Also takes a "negative input list" for denominators (Frac).
# Note: we filter n = 0 results.
def _fold_make_dca(input_list_positive: Iterable[Function],
                   input_list_negative: Iterable[Function],
                   cardinal_type: type[FunctionWithCardinality]
                   ) -> dict[float, list[Function]]:

    acd: dict[Function, float] = {} # args_cardinality_dict
    for k, lf in [(1, input_list_positive), (-1, input_list_negative)]:
        for f in lf:
            if isinstance(f, cardinal_type):
                acd.setdefault(f.f, 0)
                acd[f.f] += f.n.number * k
            else:
                acd.setdefault(f, 0)
                acd[f] += k

    dca: dict[float, list[Function]] = {} # inverse map of acd, for packing
    for f, n in acd.items():
        if n != 0:
            dca.setdefault(n, []).append(f)
    return dca

def _fold(input_list: Iterable[Function],
          variadic_type: type[FunctionVariadic],
          cardinal_type: type[FunctionWithCardinality]
    ) -> list[Function]:

    dca = _fold_make_dca(input_list, [], cardinal_type)
    ret: list[Function] = []
    for n, lf in dca.items():
        if n == 0.0:
            continue
        elif len(lf) == 1:
            f = lf[0]
            ret.append(f if n == 1.0 else cardinal_type.factory(f, _k(n)))
        elif n == 1.0:
            ret.extend(lf)
        else:
            ret.append(cardinal_type.factory(variadic_type.factory(*lf), _k(n)))
    return ret

def _k(x: float) -> Const:
    return Const.factory(x)

def _fix_const(x: Const | float) -> Const:
    return cast(Const, _fix_const_if_float(x))

def _fix_const_if_float(x: Function | float) -> Function:
    if isinstance(x, int | float):
        return _k(x)
    else:
        return x

def _fix_consts_in_iterable(it: Iterable[Function | float]) -> list[Function]:
    ret: list[Function] = []
    for x in it:
        ret.append(_fix_const_if_float(x))
    return ret

def _parens_if_sum(x: Function) -> str:
    if isinstance(x, Sum):
        return f'({x})'
    else:
        return str(x)

def _parens_if_sum_prod_frac_times(x: Function) -> str:
    if isinstance(x, Sum | Prod | Frac | Times):
        return f'({x})'
    else:
        return str(x)

_SORT_ORDER = (
    Const,
    ConstName,
    Var,
    Sum,
    Minus,
    Times,
    Prod,
    Inv,
    Frac,
    Sq,
    Sqrt,
    Pow,
    Sin,
    Cos,
    Arctan,
    ExplicitFunction,
    OpaqueFunction,
    )

@dataclass
class _DataClass(ABC):
    @abstractmethod
    def __init__(self) -> None:
        pass

# default_factory: classes should not have [] as default list value

@dataclass
class _DataSum(_DataClass):
    consts:    list[Const]    = dataclasses.field(default_factory=list)
    functions: list[Function] = dataclasses.field(default_factory=list)

@dataclass
class _DataFrac(_DataClass):
    consts: list[Const]    = dataclasses.field(default_factory=list)
    nums:   list[Function] = dataclasses.field(default_factory=list)
    dens:   list[Function] = dataclasses.field(default_factory=list)

def _merge_datas(*args: _DataClass) -> _DataClass:
    match len(args):
        case 0:
            raise ValueError('At least one argument must be passed')
        case 1:
            return args[0]

    t = type(args[0])
    for x in args[1:]:
        if type(x) is not type(args[0]):
            raise ValueError('All values must be of same type')

    ret = t()
    for field in dataclasses.fields(t):
        fn = field.name
        for x in args:
            ret.__dict__[fn].extend(x.__dict__[fn])
    return ret

def _dataclass_distribute(
        type_to_distribute: type[FunctionWithCardinality],
        n: Const | float,
        d: _DataClass,
        ) -> _DataClass:
    t = type(d)
    ret = t()
    for field in dataclasses.fields(t):
        fn = field.name
        ret.__dict__[fn].extend(
                [type_to_distribute.factory(f, _fix_const(n))
                    for f in d.__dict__[fn]])
    return ret

def _merge_datafracs_inverse(d1: _DataFrac, d2: _DataFrac) -> _DataFrac:
    ret = _DataFrac()
    ret.consts.extend(d1.consts)
    ret.consts.extend([_k(1.0/x.number) for x in d2.consts])
    ret.nums.extend(d1.nums)
    ret.nums.extend(d2.dens)
    ret.dens.extend(d1.dens)
    ret.dens.extend(d2.nums)
    return ret

def _flatten_args(t: type[Function], *args: Function) -> _DataClass:
    return _merge_datas(*[_flatten(t, f) for f in args])

def _flatten_sum(*args: Function) -> _DataSum:
    return cast(_DataSum, _flatten_args(Sum, *args))

def _flatten_prod(*args: Function) -> _DataFrac:
    return cast(_DataFrac, _flatten_args(Prod, *args))

def _flatten(t: type[Function], f: Function) -> _DataClass:

    def _flatten_variadic(t: type[Function], f: Function) -> _DataClass:
        if not isinstance(f, FunctionVariadic): raise TypeError()
        return _flatten_args(t, *f.args)

    def _flatten_default_sum(t: type[Function], f: Function) -> _DataClass:
        if t is not Sum: raise ValueError()
        return _DataSum(functions = [f])

    def _flatten_default_prod(t: type[Function], f: Function) -> _DataClass:
        if t is not Prod: raise ValueError()
        return _DataFrac(nums = [f])

    def _flatten_sum_const(t: type[Function], f: Function) -> _DataClass:
        if t is not Sum: raise ValueError()
        if not isinstance(f, Const): raise TypeError()
        return _DataSum(consts = [f])

    def _flatten_prod_const(t: type[Function], f: Function) -> _DataClass:
        if t is not Prod: raise ValueError()
        if not isinstance(f, Const): raise TypeError()
        return _DataFrac(consts = [f])

    def _flatten_sum_times(t: type[Function], f: Function) -> _DataClass:
        if t is not Sum: raise ValueError()
        if not isinstance(f, Times): raise TypeError()
        return _dataclass_distribute(Times, f.n, _flatten(t, f.f))

    def _flatten_prod_times(t: type[Function], f: Function) -> _DataClass:
        if t is not Prod: raise ValueError()
        if not isinstance(f, Times): raise TypeError()
        return _merge_datas(_DataFrac(consts = [f.n]), _flatten(t, f.f))

    def _flatten_prod_pow(t: type[Function], f: Function) -> _DataClass:
        if t is not Prod: raise ValueError()
        if not isinstance(f, Pow): raise TypeError()
        if f.n < ZERO:
            return _flatten_prod_frac(t, Frac(ONE, Pow(f.f, -f.n)))
        return _dataclass_distribute(Pow, f.n, _flatten(t, f.f))

    def _flatten_prod_frac(t: type[Function], f: Function) -> _DataClass:
        if t is not Prod: raise ValueError()
        if not isinstance(f, Frac): raise TypeError()
        return _merge_datafracs_inverse(
                cast(_DataFrac, _flatten(t, f.num)),
                cast(_DataFrac, _flatten(t, f.den)))

    d: dict[
            tuple[type[Function], type[Function]],
            Callable[[type[Function], Function], _DataClass]
        ] = {
            (Sum, Sum):     _flatten_variadic,
            (Prod, Prod):   _flatten_variadic,
            (Sum, Const):   _flatten_sum_const,
            (Prod, Const):  _flatten_prod_const,
            (Prod, Times):  _flatten_prod_times,
            (Sum, Times):   _flatten_sum_times,
            (Prod, Pow):    _flatten_prod_pow,
            (Prod, Frac):   _flatten_prod_frac,
            (Sum, Function):    _flatten_default_sum,
            (Prod, Function):   _flatten_default_prod,
        }

    for (t1, t2), func in d.items():
        if t is t1 and isinstance(f, t2):
            return func(t, f)
    raise NotImplementedError(f"flatten() has no match for {(t,type(f))}")

@cache
def _extract_any[T: Hashable](o: Any, t: type[T]) -> set[T]:
    if isinstance(o, _ExtractableMixin):
        return o.extract(t)

    s = set()
    if isinstance(o, t):
        s.add(o)
    if isinstance(o, str | bytes | bytearray):
        pass
    elif isinstance(o, list | set | tuple):
        s.update({ooo for oo in o for ooo in _extract_any(oo, t)})
    elif isinstance(o, dict):
        for collection in o.keys(), o.values():
            s.update({ooo for oo in collection for ooo in _extract_any(oo, t)})
    elif isinstance(o, Iterable):
        raise NotImplementedError('Not ready for {type(o)}')
    return s
