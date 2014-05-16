import unittest


class Test(unittest.TestCase):

    def setUp(self):
        pass

    def test_1(self):
        import slycot
        slycot.sb02mt(1,1,1,1)

    def test_2(self):
        import slycot
        from scipy import matrix
        a = matrix("-2 0.5;-1.6 -5")
        Ar, Vr, Yr, VALRr, VALDr = slycot.mb05md(a, 0.1)
