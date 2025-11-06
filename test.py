import chain as chain
from chain import *

from unittest import TestCase, main as unittest_main
import inspect

k = chain._k

class ChainTestClass(TestCase):

    def test_sort_priority(self):

        # make sure no subclass of Function is left behind in ordering
        sp = chain._SORT_ORDER
        for c in [x[1] for x in inspect.getmembers(chain, inspect.isclass)]:
            if issubclass(c, Function) and not inspect.isabstract(c):
                self.assertTrue(c in sp,
                    f'Class {c.__name__} is not in sort priority list')

        # test comparison methods

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

        self.assertFalse(X <  X)
        self.assertTrue (X <= X)
        self.assertFalse(X >  X)
        self.assertTrue (X >= X)
        self.assertTrue (X == X)
        self.assertFalse(X != X)

        self.assertTrue (ONE <  X)
        self.assertTrue (ONE <= X)
        self.assertFalse(ONE >  X)
        self.assertFalse(ONE >= X)
        self.assertFalse(ONE == X)
        self.assertTrue (ONE != X)

        self.assertFalse(X <  ONE)
        self.assertFalse(X <= ONE)
        self.assertTrue (X >  ONE)
        self.assertTrue (X >= ONE)
        self.assertFalse(X == ONE)
        self.assertTrue (X != ONE)

    def test_numerical_operators(self):
        self.assertEqual(ONE + ONE, TWO)
        self.assertEqual(X + X, Times(X, TWO))
        self.assertEqual(X + Y, Sum(X, Y))

        self.assertEqual(ONE - ONE, ZERO)
        self.assertEqual(X - X, ZERO)
        self.assertEqual(X - Y, Sum(X, Minus(Y)))

        self.assertEqual(ONE * ONE, ONE)
        self.assertEqual(X * X, Sq(X))
        self.assertEqual(X * Y, Prod(X, Y))

        self.assertEqual(ONE / ONE, ONE)
        self.assertEqual(X / X, ONE)
        self.assertEqual(X / Y, Frac(X, Y))

    def test_const(self):
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
        x = Var('x')
        self.assertTrue(X == x)
        self.assertTrue(X < Y)
        self.assertFalse(X < x)
        self.assertFalse(X == Y)
        self.assertEqual(X.der(X), ONE)
        self.assertEqual(X.der(Y), ZERO)
        self.assertEqual(str(x), 'x')

    def test_fix_consts(self):
        f1 = chain._fix_const
        f2 = chain._fix_const_if_float
        f3 = chain._fix_consts_in_iterable
        self.assertEqual(f1(1), ONE)
        self.assertEqual(f1(ONE), ONE)
        self.assertEqual(f2(1), ONE)
        self.assertEqual(f2(ONE), ONE)
        self.assertEqual(f2(X), X)
        self.assertEqual(f3(tuple()), [])
        self.assertEqual(f3((1, ONE, X)), [ONE, ONE, X])

    def test_minus(self):
        self.assertEqual(Minus.factory(X), Times.factory(X, MINUS_ONE))
        self.assertEqual(Minus.factory(Minus.factory(X)), X)

    def test_pow(self):
        three = Const.factory(3)
        self.assertEqual(
                Pow.factory(Times.factory(X, three), three),
                Times.factory(Pow.factory(X, three), Const.factory(27)))

    def test_sum(self):
        sf = Sum.factory
        self.assertEqual(sf(), ZERO)
        tests: list[tuple[list[Const | int], Const]] = [
            ([1, 2, 3], k(6)),
            ([ONE, ONE], TWO),
            ([ONE, MINUS_ONE], ZERO),
            ([k(1000), k(-1000.0)], ZERO),
            ([k(1000.1), k(-1000.0)], k(1000.1 - 1000)),
        ]
        for t in tests:
            self.assertEqual(sf(*t[0]), t[1])

        # flattening
        self.assertEqual(sf(sf(sf(sf(sf(sf(sf(1))))))), ONE)
        self.assertEqual(sf(sf(1), sf(X), sf(Y)), Sum(ONE, X, Y))
        self.assertEqual(sf(sf(1, X), sf(Y, Z)), Sum(ONE, X, Y, Z))
        msf = Minus.factory(sf(X, Y))
        self.assertEqual(sf(X, Y, Z, Minus.factory(sf(Y, Z))), X)

        # const reducing
        self.assertEqual(sf(1.0, 10, -5.5, k(8), k(-2.0)),
                         k(1.0 + 10 -5.5 + 8 -2.0))

        return
        # folding
        self.assertEqual(sf(
                X, X, X,
                Y, Minus.factory(Y),
                Z, Z,
                Times.factory(Z, 3), Times.factory(Z, -10.5),
                Minus.factory(Z), Minus.factory(Z)),
            Sum(Times(X, 3), Times(Z, 1+1+3-10.5-1-1)))
        self.assertEqual(
                sf(
                    Prod.factory(X, Y),
                    Prod.factory(Y, X)),
                Times(Prod(X, Y), 2))

        # ordering
        self.assertEqual(sf(1, X, Y), sf(Y, X, 1))

        # derivative
        self.assertEqual(
            sf(1, X, Sq(X), Y, Sq(Y)).der(X),
            sf(1, Times.factory(X, 2)))

    def test_prod_consts(self):
        tests: list[tuple[list[Const | int], Const]] = [
            ([1, 2, 3, 4], k(24)),
            ([ONE, ONE], ONE),
            ([ONE, MINUS_ONE], MINUS_ONE),
            ([ONE, 0, MINUS_ONE], ZERO),
            ([k(1000), k(-1000.0)], k(-1e6)),
            ([k(1000.1), k(-1000.0)], k(1000.1 * -1000)),
        ]
        for t in tests:
            self.assertEqual(Prod.factory(*t[0]), t[1])

    def test_prod_order(self):
        self.assertEqual(
                Prod.factory(ONE, X, Y, Z),
                Prod.factory(Z, Y, X, ONE))

    def test_prod_mult(self):
        self.assertEqual(
                Prod.factory(X, X),
                Sq.factory(X))
        self.assertEqual(
                Prod.factory(
                    Prod.factory(2, X, Y),
                    Prod.factory(2, Y, X)),
                Times.factory(Prod.factory(Sq(X), Sq(Y)), 4))
        self.assertEqual(
                Prod.factory(
                    X,
                    Pow.factory(X, -2)),
                Inv(X))

    def test_prod_frac(self):
        pf = Prod.factory
        a = Var('a')
        b = Var('b')
        c = Var('c')
        d = Var('d')
        e = Var('e')
        f = Var('f')
        f1 = Frac(X, Y)
        f2 = Frac(a, b)
        f3 = Frac(c, d)
        f4 = Frac(e, f)
        self.assertEqual(pf(f1, f2, f3, f4), Frac(Prod(a,c,e,X),Prod(b,d,f,Y)))
        self.assertEqual(pf(Frac(f1, f2), Frac(f3, f4)),
                         Frac(Prod(b,c,f,X),Prod(a,d,e,Y)))
        self.assertEqual(pf(Sq(X), Frac(ONE,X), Frac(ONE, X)), ONE)
        self.assertEqual(pf(Y, Frac(ONE,Inv(X))), Prod(X, Y))
        self.assertEqual(pf(Y, Frac(ONE, X), Frac(ONE,Inv(X))), Y)
        self.assertEqual(pf(Inv(Y), Frac(ONE,X)), Inv(Prod(X, Y)))

    def test_prod_der(self):
        a = X
        b = Sum.factory(Times(X, TWO), Y)
        c = Z
        self.assertEqual(
                Prod.factory(a, b, c).der(X),
                Sum.factory(
                    Prod.factory(b, c),
                    Prod.factory(a, c, TWO),
                    ))

    def test_frac_factory(self):
        ff = Frac.factory

        self.assertEqual(ff(25, -5), -5)
        self.assertEqual(ff(ZERO, ONE), ZERO)
        self.assertRaises(ZeroDivisionError, ff, X, ZERO)
        self.assertEqual(ff(X, ONE), X)
        self.assertEqual(ff(ONE, TWO), k(.5))
        self.assertEqual(ff(ONE, X), Inv(X))
        self.assertEqual(ff(TWO, X), Times(Inv(X), TWO))
        self.assertEqual(ff(X, TWO), Times(X, k(.5)))

        pz = Prod(X, Prod(Y, ZERO, Z), Sum(X, Y))
        self.assertEqual(ff(pz, X), ZERO)
        self.assertRaises(ZeroDivisionError, ff, X, pz)

        self.assertEqual(ff(X, X), ONE)
        self.assertEqual(ff(Times.factory(X, 2), X), TWO)
        self.assertEqual(ff(Sq(X), X), X)
        self.assertEqual(ff(X, Sq(X)), Inv(X))
        self.assertEqual(ff(Sq(X), Sq(Y)), Sq(Frac(X, Y)))

        self.assertEqual(ff(X, Minus(X)), MINUS_ONE)
        self.assertEqual(ff(Times.factory(X, -2), X), k(-2))
        self.assertEqual(ff(Sq(X), Minus(X)), Minus(X))
        self.assertEqual(ff(Minus(X), Sq(X)), Minus(Inv(X)))
        self.assertEqual(ff(Minus(Sq(X)), Sq(Y)), Minus(Sq(Frac(X, Y))))

        a = Var('a')
        b = Var('b')
        c = Var('c')
        d = Var('d')
        e = Var('e')
        f = Var('f')
        f1 = Frac(X, Y)
        f2 = Frac(a, b)
        f3 = Frac(c, d)
        f4 = Frac(e, f)
        self.assertEqual(ff(Frac(f1, f2), Frac(f3, f4)),
                         Frac(Prod(b, d, e, X), Prod(a, c, f, Y)))
        self.assertEqual(ff(Sq(Frac(f1, f2)), Frac(f3, f4)),
                         Prod.factory(
                             Sq(Frac(f1, f2)),
                             Frac(Prod(d, e), Prod(c, f))))
        self.assertEqual(ff(Frac(Sq(f1), Sq(f2)), Frac(f3, f4)),
                         Prod.factory(
                             Sq(Frac(f1, f2)),
                             Frac(Prod(d, e), Prod(c, f))))
        self.assertEqual(ff(  Sum(a, b, c), Sum(a, b, d)),
                         Frac(Sum(a, b, c), Sum(a, b, d)))
        self.assertEqual(ff(
            Prod(Times(Pow(X, k(1.5)), k(2.5)), Pow(Frac(X, Y), k(-4))),
            X),
            Times(
                Frac(
                    Pow(Var('y'), Const(4.0)),
                    Pow(Var('x'), Const(3.5))),
                Const(2.5)))
        self.assertEqual(ff(Sq(X), Prod(Sq(Y), Z)),
                         Prod.factory(Sq(Frac(X,Y)), Inv(Z)))
        self.assertEqual(ff(Sq(X), Prod(Sq(Y), Pow(Z, k(3)))),
                         Prod.factory(Sq(Frac(X,Y)), Pow(Z, k(-3))))

    def test_frac_der(self):
        num = Sum.factory(Sq(X), Y, 1)
        den = Sum.factory(Times(Sq(X), k(3)), X, Sq(Y), Z)
        self.assertEqual(
            Frac.factory(num, den).der(X),
                Frac.factory(
                    Sum.factory(
                        Prod.factory(Times(X, TWO), den),
                        Minus.factory(
                            Prod.factory(num, Sum.factory(Times(X, k(6)), 1)),
                    )),
                    Sq.factory(den))
        )

    def test_frac_str(self):
        sf = Sum.factory
        pf = Prod.factory
        ff = Frac.factory
        self.assertEqual(str(Frac(ONE, ONE)), '1/1')
        self.assertEqual(str(Frac(Sum(ONE, ONE), ONE)), '(1 + 1)/1')
        self.assertEqual(str(Frac(Prod(ONE, ONE), ONE)), '1*1/1')
        self.assertEqual(str(Frac(ONE, Sum(ONE, ONE))), '1/(1 + 1)')
        self.assertEqual(str(Frac(ONE, Prod(ONE, ONE))), '1/(1*1)')

    def test_sin_cos(self):
        sf = Sin.factory
        sx = sf(X)
        cf = Cos.factory
        cx = cf(X)
        self.assertEqual(sx.der(X), cx)
        self.assertEqual(cx.der(X), Minus(sx))
        self.assertEqual(sf(cx).der(X), Minus(Prod(sx, cf(cx))))
        self.assertEqual(cf(sx).der(X), Minus(Prod(sf(sx), cx)))

        self.assertEqual(sf(ZERO), ZERO)
        self.assertEqual(sf(PI), ZERO)
        pi = math.pi
        lst = [
                (pi/2, ONE),
                (3*pi/2, MINUS_ONE),
                (5*pi/2, ONE),
                (7*pi/2, MINUS_ONE),
                (-pi/2, MINUS_ONE),
                (-3*pi/2, ONE),
                (-pi, ZERO),
                (-2*pi, ZERO),
                (2*pi, ZERO),
                (3*pi, ZERO),
        ]
        for x, y in lst:
            self.assertEqual(sf(k(x)), y)

        self.assertEqual(cf(ZERO), ONE)
        self.assertEqual(cf(PI), MINUS_ONE)
        pi = math.pi
        lst = [
                (pi/2, ZERO),
                (3*pi/2, ZERO),
                (5*pi/2, ZERO),
                (7*pi/2, ZERO),
                (-pi/2, ZERO),
                (-3*pi/2, ZERO),
                (-pi, MINUS_ONE),
                (-2*pi, ONE),
                (2*pi, ONE),
                (3*pi, MINUS_ONE),
        ]
        for x, y in lst:
            self.assertEqual(cf(k(x)), y)

    def test_explicit_function(self):
        name = 'cosxy'
        f = Cos.factory(Sum(Times(X,k(2)),Times(Y,k(3))))
        uf = ExplicitFunction(name, f)
        dufx = ExplicitFunction(name, f.der(X), (X,))
        dufy = ExplicitFunction(name, f.der(Y), (Y,))
        dufxy = ExplicitFunction(name, f.der(X).der(Y), (X,Y))
        dufyx = ExplicitFunction(name, f.der(Y).der(X), (Y,X))
        self.assertEqual(uf.der(X), dufx)
        self.assertEqual(uf.der(Y), dufy)
        self.assertEqual(uf.der(X).der(Y), dufxy)
        self.assertEqual(uf.der(Y).der(X), dufyx)
        self.assertEqual(dufx.f, f.der(X))
        self.assertEqual(dufy.f, f.der(Y))
        self.assertEqual(dufxy.f, f.der(X).der(Y))
        self.assertEqual(dufyx.f, f.der(Y).der(X))
        self.assertEqual(uf.der(Z), ZERO)
        self.assertEqual(dufx.der(Z), ZERO)
        self.assertEqual(dufy.der(Z), ZERO)
        self.assertEqual(dufxy.der(Z), ZERO)
        self.assertEqual(dufyx.der(Z), ZERO)

    def test_opaque_function(self):
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
        self.assertEqual(f.der(Z), ZERO)
        self.assertEqual(dfx.der(Z), ZERO)
        self.assertEqual(dfy.der(Z), ZERO)
        self.assertEqual(dfxy.der(Z), ZERO)
        self.assertEqual(dfyx.der(Z), ZERO)

    def test_arctan(self):
        self.assertEqual(Arctan.factory(X).der(X), Inv(Sum.factory(ONE, Sq(X))))

if __name__ == '__main__':
    unittest_main()
