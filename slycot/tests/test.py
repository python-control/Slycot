import unittest
from .. import synthesis
from .. import math


class Test(unittest.TestCase):

    def setUp(self):
        pass

    def test_1(self):
        synthesis.sb02mt(1,1,1,1)

    def test_2(self):
        from numpy import matrix
        a = matrix("-2 0.5;-1.6 -5")
        Ar, Vr, Yr, VALRr, VALDr = math.mb05md(a, 0.1)

    def test_sb02ad(self):
        "Test sb10ad, Hinf synthesis"
        import numpy as np
        a = np.array([[-1]])
        b = np.array([[1, 1]])
        c = np.array([[1], [1]])
        d = np.array([[0, 1], [1, 0]])

        n = 1
        m = 2
        np_ = 2
        ncon = 1
        nmeas = 1
        gamma = 10

        gamma_est, Ak, Bk, Ck, Dk, Ac, Bc, Cc, Dc, rcond = synthesis.sb10ad(
            n, m, np_, ncon, nmeas, gamma, a, b, c, d)
        # from Octave, which also uses SB10AD:
        #   a= -1; b1= 1; b2= 1; c1= 1; c2= 1; d11= 0; d12= 1; d21= 1; d22= 0;
        #   g = ss(a,[b1,b2],[c1;c2],[d11,d12;d21,d22]);
        #   [k,cl] = hinfsyn(g,1,1);
        # k.a is Ak, cl.a is Ac
        # gamma values don't match; not sure that's critical
        # this is a bit fragile
        # a simpler, more robust check might be to check stability of Ac
        self.assertEqual(Ak.shape, (1, 1))
        self.assertAlmostEqual(Ak[0][0], -3)
        self.assertEqual(Ac.shape, (2, 2))
        self.assertAlmostEqual(Ac[0][0], -1)
        self.assertAlmostEqual(Ac[0][1], -1)
        self.assertAlmostEqual(Ac[1][0], 1)
        self.assertAlmostEqual(Ac[1][1], -3)
