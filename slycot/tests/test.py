import unittest
from .. import synthesis
from .. import math


class Test(unittest.TestCase):

    def setUp(self):
        pass

    def test_1(self):
        synthesis.sb02mt(1,1,1,1)

    def test_2(self):
        from scipy import matrix
        a = matrix("-2 0.5;-1.6 -5")
        Ar, Vr, Yr, VALRr, VALDr = math.mb05md(a, 0.1)
