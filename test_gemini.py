#!/usr/bin/env python3

"""
Comprehensive unit test file for chain.py.
This file tests standard operations, factory simplification logic,
error handling, special derivative cases, and internal logic.
"""

# Import all classes and functions to be tested
import chain as chain
from chain import *
from chain import (
    _k, _fix_const, _fix_const_if_float, _fix_consts_in_iterable,
    _flatten_sum, _flatten_prod, _DataSum, _DataFrac, _merge_datas
)

import math
import inspect
from unittest import TestCase, main as unittest_main

# --- Factory Aliases ---
k = chain._k
sf = Sum.factory
pf = Prod.factory
tf = Times.factory
minusf = Minus.factory
fracf = Frac.factory
powf = Pow.factory
sqf = Sq.factory
sqrtf = Sqrt.factory
invf = Inv.factory
sinf = Sin.factory
cosf = Cos.factory
arctanf = Arctan.factory


# --- Test Variables ---
X = Var('x')
Y = Var('y')
Z = Var('z')
A = Var('a')
B = Var('b')
C = Var('c')

class GeminiTestClass(TestCase):

    # ======================================================================
    # == Setup and Basic Class Tests (from test.py)
    # ======================================================================

    def test_sort_priority(self):
        """
        Tests that all Function subclasses are in the _SORT_ORDER list
        and that comparison magic methods work correctly.
        """
        sp = chain._SORT_ORDER
        all_funcs = [x[1] for x in inspect.getmembers(chain, inspect.isclass)]

        for c_ in all_funcs:
            if issubclass(c_, Function) and not inspect.isabstract(c_):
                self.assertTrue(c_ in sp,
                    f'Class {c_.__name__} is not in sort priority list')

        # Test comparison methods
        self.assertLess(sp.index(Const), sp.index(Var))
        self.assertTrue (X <  Y)
        self.assertTrue (X <= Y)
        self.assertFalse(X >  Y)
        self.assertFalse(X >= Y)
        self.assertFalse(X == Y)
        self.assertTrue (X != Y)

        self.assertFalse(Y <  X)
        self.assertFalse(Y <= X)
        self.assertTrue (Y >  X)
        self.assertTrue (Y >= X)
        self.assertFalse(Y == X)
        self.assertTrue (Y != X)

        self.assertTrue (X == X)
        self.assertFalse(X != X)
        self.assertTrue (ONE <  X)
        self.assertFalse(X <  ONE)

    def test_numerical_operators(self):
        """Tests the overloaded numerical operators (+, -, *, /)."""
        self.assertEqual(ONE + ONE, TWO)
        self.assertEqual(X + X, tf(X, TWO))
        self.assertEqual(X + Y, sf(X, Y))

        self.assertEqual(ONE - ONE, ZERO)
        self.assertEqual(X - X, ZERO)
        self.assertEqual(X - Y, sf(X, minusf(Y)))

        self.assertEqual(ONE * ONE, ONE)
        self.assertEqual(X * X, sqf(X))
        self.assertEqual(X * Y, pf(X, Y))

        self.assertEqual(ONE / ONE, ONE)
        self.assertEqual(X / X, ONE)
        self.assertEqual(X / Y, fracf(X, Y))

    def test_const(self):
        """Tests the Const class and its factory."""
        self.assertEqual(k(1), k(1.0))
        self.assertEqual(ONE.der(X), ZERO)
        self.assertEqual(k(True), ONE)
        self.assertEqual(k(False), ZERO)
        self.assertTrue(ONE == Const(1))
        self.assertTrue(ONE < TWO)
        self.assertFalse(ONE < Const(1))
        self.assertFalse(ONE == TWO)
        self.assertEqual(str(ONE), str(1))
        self.assertEqual(str(k(1.5)), str(1.5))
        for x in range(-5, 5):
            self.assertEqual(k(x).der(X), ZERO)

    def test_variable(self):
        """Tests the Var class."""
        x = Var('x')
        self.assertTrue(X == x)
        self.assertTrue(X < Y)
        self.assertFalse(X < x)
        self.assertFalse(X == Y)
        self.assertEqual(X.der(X), ONE)
        self.assertEqual(X.der(Y), ZERO)
        self.assertEqual(str(x), 'x')

    def test_fix_consts(self):
        """Tests internal _fix_const helper functions."""
        f1 = _fix_const
        f2 = _fix_const_if_float
        f3 = _fix_consts_in_iterable
        self.assertEqual(f1(1), ONE)
        self.assertEqual(f1(ONE), ONE)
        self.assertEqual(f2(1), ONE)
        self.assertEqual(f2(ONE), ONE)
        self.assertEqual(f2(X), X)
        self.assertEqual(f3(tuple()), [])
        self.assertEqual(f3((1, ONE, X)), [ONE, ONE, X])

    # ======================================================================
    # == Systematic Error Handling Tests (High Priority)
    # ======================================================================

    def test_error_handling(self):
        """Tests that factories raise exceptions on invalid inputs."""
        # 0^0
        self.assertRaises(ValueError, powf, 0, 0)
        # 0^(-1)
        self.assertRaises(ValueError, powf, 0, -1)
        # (-4)^0.5 (complex result)
        self.assertRaises(ValueError, powf, k(-4), 0.5)

        # 1/0
        self.assertRaises(ValueError, invf, 0)
        self.assertRaises(ValueError, invf, ZERO)

        # X/0
        self.assertRaises(ZeroDivisionError, fracf, X, 0)
        self.assertRaises(ZeroDivisionError, fracf, X, ZERO)

        # Invalid cardinality
        self.assertRaises(ValueError, minusf, X, 2)
        self.assertRaises(ValueError, sqf, X, 3)
        self.assertRaises(ValueError, invf, X, 2)
        self.assertRaises(ValueError, sqrtf, X, 2)

    # ======================================================================
    # == __repr__ and __hash__ Tests (Low Priority)
    # ======================================================================

    def test_repr_and_hash(self):
        """Tests __repr__ and __hash__ for consistency."""
        self.assertEqual(repr(X), "Var('x')")
        self.assertEqual(repr(ONE), "Const(1.0)")
        self.assertEqual(repr(sf(X, 1)), "Sum(Const(1.0), Var('x'))")

        # Test hash consistency
        x_prime = Var('x')
        one_prime = k(1)
        self.assertEqual(hash(X), hash(x_prime))
        self.assertEqual(hash(ONE), hash(one_prime))
        self.assertNotEqual(hash(X), hash(Y))

        # Test usage in a dictionary (relies on hash and eq)
        d = {X: 'foo', ONE: 'bar'}
        self.assertEqual(d[x_prime], 'foo')
        self.assertEqual(d[one_prime], 'bar')
        self.assertTrue(Y not in d)

    # ======================================================================
    # == Internal Flatten Logic Tests (Medium Priority)
    # ======================================================================

    def test_internal_flatten_sum(self):
        """Tests the internal _flatten_sum logic."""
        # _flatten_sum(sf(X, tf(sf(Y, 1), 2)))
        # -> _flatten_sum(X, tf(Y, 2), tf(1, 2))
        # -> DataSum(consts=[k(2)], functions=[X, tf(Y, 2)])
        d = _flatten_sum(sf(X, tf(sf(Y, 1), 2)))
        self.assertEqual(d.consts, [k(2)])
        self.assertEqual(set(d.functions), {X, tf(Y, 2)})

    def test_internal_flatten_prod(self):
        """Tests the internal _flatten_prod logic."""
        # _flatten_prod(pf(X, powf(tf(Y, 2), 3), fracf(Z, 4)))
        # -> X, (Y*2)^3, Z/4
        # -> X, (Y^3 * 2^3), Z * 1/4
        # -> X, Y^3, 8, Z, 0.25
        # -> DataFrac(consts=[8, 0.25], nums=[X, powf(Y, 3), Z], dens=[])
        d = _flatten_prod(pf(X, powf(tf(Y, 2), 3), fracf(Z, 4)))
        self.assertEqual(set(d.consts), {k(2)})
        self.assertEqual(set(d.nums), {X, powf(Y, 3), Z})
        self.assertEqual(d.dens, [])

        # Test Frac flattening
        # _flatten_prod(fracf(X, pf(Y, Z)))
        # -> DataFrac(consts=[], nums=[X], dens=[Y, Z])
        d2 = _flatten_prod(fracf(X, pf(Y, Z)))
        self.assertEqual(d2.consts, [])
        self.assertEqual(d2.nums, [X])
        self.assertEqual(set(d2.dens), {Y, Z})

    # ======================================================================
    # == Sum / Minus Factory Tests
    # ======================================================================

    def test_minus_factory(self):
        """Tests Minus.factory simplifications."""
        self.assertEqual(minusf(X), tf(X, MINUS_ONE))
        self.assertEqual(minusf(minusf(X)), X)
        self.assertEqual(minusf(tf(X, 2)), tf(X, -2))
        self.assertEqual(minusf(k(5)), k(-5))

    def test_sum_factory(self):
        """Tests Sum.factory logic (const folding, flattening, folding)."""
        # Empty
        self.assertEqual(sf(), ZERO)

        # Const folding
        tests: list[tuple[list[Const | int], Const]] = [
            ([1, 2, 3], k(6)),
            ([ONE, ONE], TWO),
            ([ONE, MINUS_ONE], ZERO),
            ([k(1000), k(-1000.0)], ZERO),
            ([k(1000.1), k(-1000.0)], k(1000.1 - 1000)),
        ]
        for t in tests:
            self.assertEqual(sf(*t[0]), t[1])

        # Flattening
        self.assertEqual(sf(sf(sf(sf(sf(sf(sf(1))))))), ONE)
        self.assertEqual(sf(sf(1), sf(X), sf(Y)), sf(ONE, X, Y))
        self.assertEqual(sf(sf(1, X), sf(Y, Z)), sf(ONE, X, Y, Z))
        msf = minusf(sf(X, Y))
        self.assertEqual(sf(X, Y, Z, minusf(sf(Y, Z))), X)
        self.assertEqual(sf(X, 1, Y, -1, Z), sf(X, Y, Z))

        # Folding (this was disabled in test.py)
        self.assertEqual(sf(X, X), tf(X, 2))
        self.assertEqual(sf(X, X, X), tf(X, 3))
        self.assertEqual(sf(X, minusf(X)), ZERO)
        self.assertEqual(sf(X, X, Y, minusf(Y)), tf(X, 2))
        self.assertEqual(sf(pf(X,Y), pf(Y,X)), tf(pf(X,Y), 2))

        self.assertEqual(sf(
                X, X, X,  # 3*X
                Y, minusf(Y), # 0
                Z, Z, # 2*Z
                tf(Z, 3), # 3*Z
                tf(Z, -10.5), # -10.5*Z
                minusf(Z), minusf(Z)), # -2*Z
            sf(tf(X, 3), tf(Z, 2 + 3 - 10.5 - 2)))

        # Ordering
        self.assertEqual(sf(1, X, Y), sf(Y, X, 1))

        # Derivative
        self.assertEqual(
            sf(1, X, sqf(X), Y, sqf(Y)).der(X),
            sf(1, tf(X, 2)))

    # ======================================================================
    # == Prod / Times Factory Tests
    # ======================================================================

    def test_prod_factory_consts(self):
        """Tests Prod.factory with only constants."""
        tests: list[tuple[list[Const | int], Const]] = [
            ([1, 2, 3, 4], k(24)),
            ([ONE, ONE], ONE),
            ([ONE, MINUS_ONE], MINUS_ONE),
            ([ONE, 0, MINUS_ONE], ZERO),
            ([k(1000), k(-1000.0)], k(-1e6)),
            ([k(1000.1), k(-1000.0)], k(1000.1 * -1000)),
        ]
        for t in tests:
            self.assertEqual(pf(*t[0]), t[1])

    def test_prod_factory_simplifications(self):
        """Tests Prod.factory logic (folding, flattening)."""
        # Empty
        self.assertEqual(pf(), ONE)

        # Zero
        self.assertEqual(pf(X, Y, Z, 0), ZERO)

        # Flattening
        self.assertEqual(pf(pf(X, Y), pf(Z, A)), pf(A, X, Y, Z))
        self.assertEqual(pf(pf(X, 2), pf(Y, 3)), tf(pf(X, Y), 6))

        # Folding (Powers)
        self.assertEqual(pf(X, X), sqf(X))
        self.assertEqual(pf(X, X, X), powf(X, 3))
        self.assertEqual(pf(X, invf(X)), ONE)
        self.assertEqual(pf(X, X, invf(X)), X)
        self.assertEqual(pf(X, powf(X, 2)), powf(X, 3))
        self.assertEqual(pf(X, powf(X, -2)), invf(X))

        self.assertEqual(
                pf(
                    pf(2, X, Y),
                    pf(2, Y, X)),
                tf(pf(sqf(X), sqf(Y)), 4))

        # Ordering
        self.assertEqual(pf(ONE, X, Y, Z), pf(Z, Y, X, ONE))

        # Frac/Inv folding
        self.assertEqual(pf(fracf(X, Y), fracf(Y, X)), ONE)
        self.assertEqual(pf(fracf(X, Y), fracf(A, B)), fracf(pf(A, X), pf(B, Y)))
        self.assertEqual(pf(fracf(X, Y), Y), X)
        self.assertEqual(pf(fracf(X, Y), invf(X)), invf(Y))
        self.assertEqual(pf(Y, fracf(ONE,invf(X))), pf(X, Y))
        self.assertEqual(pf(Y, fracf(ONE, X), fracf(ONE,invf(X))), Y)
        self.assertEqual(pf(invf(Y), fracf(ONE,X)), invf(pf(X, Y)))

    def test_prod_derivative(self):
        """Tests the product rule derivative."""
        a = X
        b = sf(tf(X, TWO), Y) # 2X + Y
        c = Z

        # d(a*b*c)/dx = (da/dx)*b*c + a*(db/dx)*c + a*b*(dc/dx)
        #             = 1 * (2X + Y) * Z + X * 2 * Z + X * (2X + Y) * 0
        self.assertEqual(
                pf(a, b, c).der(X),
                sf(
                    pf(b, c),     # (2X + Y) * Z
                    pf(a, c, TWO) # X * Z * 2
                    ))

    def test_times_factory(self):
        """Tests Times.factory simplifications."""
        self.assertEqual(tf(X, 0), ZERO)
        self.assertEqual(tf(X, 1), X)
        self.assertEqual(tf(X, -1), minusf(X))
        self.assertEqual(tf(k(5), 2), k(10))

        # Flattening
        self.assertEqual(tf(tf(X, 2), 3), tf(X, 6))

        # Interaction with Minus
        self.assertEqual(tf(minusf(X), 2), tf(X, -2))
        self.assertEqual(tf(minusf(X), -2), tf(X, 2))

    # ======================================================================
    # == Pow / Sq / Sqrt / Inv Factory Tests (High Priority)
    # ======================================================================

    def test_pow_factory(self):
        """Tests Pow.factory simplifications."""
        self.assertEqual(powf(X, 1), X)
        self.assertEqual(powf(X, 0), ONE)
        self.assertEqual(powf(ONE, 5), ONE)
        self.assertEqual(powf(TWO, 3), k(8))

        # Flattening
        self.assertEqual(powf(powf(X, 2), 3), powf(X, 6))

        # Distribution
        self.assertEqual(powf(tf(X, 2), 3), tf(powf(X, 3), 8))
        self.assertEqual(powf(tf(X, 2), -1), tf(invf(X), 0.5))

        # Special types
        self.assertEqual(powf(X, 2), sqf(X))
        self.assertEqual(powf(X, 0.5), sqrtf(X))
        self.assertEqual(powf(X, -1), invf(X))

    def test_sq_factory(self):
        """Tests Sq.factory simplifications."""
        self.assertEqual(sqf(k(3)), k(9))
        self.assertEqual(sqf(powf(X, 3)), powf(X, 6))
        # Ensure negative Times is handled
        self.assertEqual(sqf(tf(X, -2)), tf(sqf(X), 4))
        self.assertEqual(sqf(minusf(X)), sqf(X))

    def test_sqrt_factory(self):
        """Tests Sqrt.factory simplifications."""
        self.assertEqual(sqrtf(k(9)), k(3))
        self.assertEqual(sqrtf(powf(X, 4)), sqf(X))
        self.assertEqual(sqrtf(sqf(X)), X) # Note: This assumes X > 0

    def test_inv_factory(self):
        """Tests Inv.factory simplifications."""
        self.assertEqual(invf(k(4)), k(0.25))
        self.assertEqual(invf(invf(X)), X)
        self.assertEqual(invf(powf(X, 2)), powf(X, -2))
        self.assertEqual(invf(powf(X, -3)), powf(X, 3))
        self.assertEqual(invf(fracf(X, Y)), fracf(Y, X))

    # ======================================================================
    # == Frac Factory Tests
    # ======================================================================

    def test_frac_factory(self):
        """Tests Frac.factory simplifications."""
        self.assertEqual(fracf(25, -5), k(-5))
        self.assertEqual(fracf(ZERO, ONE), ZERO)
        self.assertEqual(fracf(X, ONE), X)
        self.assertEqual(fracf(ONE, TWO), k(.5))
        self.assertEqual(fracf(ONE, X), invf(X))
        self.assertEqual(fracf(TWO, X), tf(invf(X), TWO))
        self.assertEqual(fracf(X, TWO), tf(X, k(.5)))

        # Identity
        self.assertEqual(fracf(X, X), ONE)

        # Power reduction
        self.assertEqual(fracf(tf(X, 2), X), TWO)
        self.assertEqual(fracf(sqf(X), X), X)
        self.assertEqual(fracf(X, sqf(X)), invf(X))

        # Minus handling
        self.assertEqual(fracf(X, minusf(X)), MINUS_ONE)
        self.assertEqual(fracf(tf(X, -2), X), k(-2))
        self.assertEqual(fracf(sqf(X), minusf(X)), minusf(X))
        self.assertEqual(fracf(minusf(X), sqf(X)), minusf(invf(X)))

        # Recursive Frac
        self.assertEqual(fracf(fracf(X, Y), fracf(A, B)), fracf(pf(B, X), pf(A, Y)))

        # Complex reduction
        self.assertEqual(fracf(sqf(X), pf(sqf(Y), Z)),
                         pf(sqf(fracf(X,Y)), invf(Z)))
        self.assertEqual(fracf(sqf(X), pf(sqf(Y), powf(Z, 3))),
                         pf(sqf(fracf(X,Y)), powf(Z, -3)))

    def test_frac_derivative(self):
        """Tests the quotient rule derivative."""
        num = sf(sqf(X), Y, 1) # X^2 + Y + 1
        den = sf(tf(X, 3), Z) # 3X + Z

        # d(num/den)/dx = (num'*den - num*den') / den^2
        # num' = 2X
        # den' = 3

        self.assertEqual(
            fracf(num, den).der(X),
                fracf(
                    sf(
                        pf(tf(X, TWO), den), # (2X) * (3X + Z)
                        minusf(
                            pf(num, k(3)), # (X^2 + Y + 1) * 3
                    )),
                    sqf(den)) # (3X + Z)^2
        )

    def test_frac_str(self):
        """Tests __str__ method for Frac to ensure parentheses."""
        self.assertEqual(str(fracf(X, Y)), 'x/y')
        self.assertEqual(str(fracf(sf(X, Y), Z)), '(x + y)/z')
        self.assertEqual(str(fracf(pf(X, Y), Z)), 'x*y/z')
        self.assertEqual(str(fracf(X, sf(Y, Z))), 'x/(y + z)')
        self.assertEqual(str(fracf(X, pf(Y, Z))), 'x/(y*z)')
        self.assertEqual(str(fracf(ONE, minusf(X))), '-1/x')

    # ======================================================================
    # == Trig / User / Opaque Function Tests
    # ======================================================================

    def test_sin_cos_factory(self):
        """Tests Sin/Cos factory constant folding."""
        # Sin
        self.assertEqual(sinf(ZERO), ZERO)
        self.assertEqual(sinf(PI), ZERO)
        pi = math.pi
        lst_sin = [
                (pi/2, ONE), (3*pi/2, MINUS_ONE),
                (-pi/2, MINUS_ONE), (-pi, ZERO), (2*pi, ZERO),
        ]
        for x, y in lst_sin:
            self.assertEqual(sinf(k(x)), y)

        # Cos
        self.assertEqual(cosf(ZERO), ONE)
        self.assertEqual(cosf(PI), MINUS_ONE)
        lst_cos = [
                (pi/2, ZERO), (3*pi/2, ZERO),
                (-pi/2, ZERO), (-pi, MINUS_ONE), (2*pi, ONE),
        ]
        for x, y in lst_cos:
            self.assertEqual(cosf(k(x)), y)

    def test_trig_derivatives(self):
        """Tests Sin/Cos derivatives."""
        sx = sinf(X)
        cx = cosf(X)
        self.assertEqual(sx.der(X), cx)
        self.assertEqual(cx.der(X), minusf(sx))

        # Chain rule
        self.assertEqual(sinf(cx).der(X), minusf(pf(sx, cosf(cx))))
        self.assertEqual(cosf(sx).der(X), minusf(pf(sinf(sx), cx)))

    def test_arctan_factory_and_der(self):
        """Tests Arctan factory and derivative."""
        self.assertEqual(arctanf(ZERO), ZERO)
        self.assertEqual(arctanf(ONE), k(math.pi/4))
        self.assertEqual(arctanf(MINUS_ONE), k(-math.pi/4))

        self.assertEqual(arctanf(X).der(X), invf(sf(ONE, sqf(X))))

        # Chain rule
        self.assertEqual(arctanf(sqf(X)).der(X),
            pf(invf(sf(ONE, powf(X, 4))), tf(X, TWO)))

    def test_user_function(self):
        """Tests UserFunction behavior."""
        name = 'cosxy'
        f = cosf(sf(tf(X,k(2)),tf(Y,k(3))))
        uf = UserFunction(name, f)

        # Check derivatives
        dufx = UserFunction(name, f.der(X), (X,))
        dufy = UserFunction(name, f.der(Y), (Y,))
        dufxy = UserFunction(name, f.der(X).der(Y), (X,Y))
        dufyx = UserFunction(name, f.der(Y).der(X), (Y,X))

        self.assertEqual(uf.der(X), dufx)
        self.assertEqual(uf.der(Y), dufy)
        self.assertEqual(uf.der(X).der(Y), dufxy)
        self.assertEqual(uf.der(Y).der(X), dufyx)

        # Check underlying function
        self.assertEqual(dufx.f, f.der(X))
        self.assertEqual(dufy.f, f.der(Y))
        self.assertEqual(dufxy.f, f.der(X).der(Y))
        self.assertEqual(dufyx.f, f.der(Y).der(X))

        # Check derivative w.r.t. unrelated var
        self.assertEqual(uf.der(Z), ZERO)
        self.assertEqual(dufx.der(Z), ZERO)
        self.assertEqual(dufyx.der(Z), ZERO)

    def test_opaque_function(self):
        """Tests OpaqueFunction behavior."""
        name = 'opaque'
        xy = (X, Y)
        f = OpaqueFunction(name, xy)

        dfx = OpaqueFunction(name, xy, (X,))
        dfy = OpaqueFunction(name, xy, (Y,))
        dfxy = OpaqueFunction(name, xy, (X,Y))
        dfyx = OpaqueFunction(name, xy, (Y,X))

        self.assertEqual(f.der(X), dfx)
        self.assertEqual(f.der(Y), dfy)
        self.assertEqual(f.der(X).der(Y), dfxy)
        self.assertEqual(f.der(Y).der(X), dfyx)

        # Check derivative w.r.t. unrelated var
        self.assertEqual(f.der(Z), ZERO)
        self.assertEqual(dfx.der(Z), ZERO)
        self.assertEqual(dfyx.der(Z), ZERO)

    # ======================================================================
    # == Derivative Special Case Tests (High Priority)
    # ======================================================================

    def test_derivative_special_cases(self):
        """
        Tests special derivative optimizations (e.g., sin^2, cos^2, sqrt).
        """

        # 1. d(sin^2(f))/dx = 2*sin(f)*cos(f)*f' = sin(2f)*f'
        f_sin = sqf(sinf(X))
        # d/dx = sin(2X) * 1
        self.assertEqual(f_sin.der(X), sinf(tf(X, 2)))

        # Test with chain rule
        f_sin_chain = sqf(sinf(sqf(X)))
        # d/dx = sin(2*X^2) * (2X)
        self.assertEqual(f_sin_chain.der(X),
            pf(sinf(tf(sqf(X), TWO)), tf(X, TWO)))

        # 2. d(cos^2(f))/dx = 2*cos(f)*(-sin(f))*f' = -sin(2f)*f'
        f_cos = sqf(cosf(X))
        # d/dx = -sin(2X) * 1
        self.assertEqual(f_cos.der(X), minusf(sinf(tf(X, 2))))

        # Test with chain rule
        f_cos_chain = sqf(cosf(sqf(X)))
        # d/dx = -sin(2*X^2) * (2X)
        self.assertEqual(f_cos_chain.der(X),
            minusf(pf(sinf(tf(sqf(X), TWO)), tf(X, TWO))))

if __name__ == '__main__':
    unittest_main()
