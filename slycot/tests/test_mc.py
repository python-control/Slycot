#!/usr/bin/env python
#
# test_mc.py - test suite for polynomial and rational function manipulation
# bnavigator <code@bnavigator.de>, Aug 2019

import unittest
import warnings

from slycot import mc01td


class test_mc(unittest.TestCase):

    def test_mc01td(self):
        """ test_mc01td: doc example
        data from http://slicot.org/objects/software/shared/doc/MC01TD.html
        """
        (dp, stable, nz) = mc01td('C', 4, [2, 0, 1, -1, 1])
        self.assertEqual(dp, 4)
        self.assertEqual(stable, 0)
        self.assertEqual(nz, 2)

    def test_mc01td_D(self):
        """ test_mc01td_D: test discrete option """
        (dp, stable, nz) = mc01td('D', 3, [1, 2, 3, 4])
        self.assertEqual(dp, 3)
        self.assertEqual(stable, 1)
        self.assertEqual(nz, 0)
        (dp, stable, nz) = mc01td('D', 3, [4, 3, 2, 1])
        self.assertEqual(dp, 3)
        self.assertEqual(stable, 0)
        self.assertEqual(nz, 3)

    def test_mc01td_warnings(self):
        """ test_mc01td_warnings: Test warnings """
        T = [([0, 0], "entry P(x) is the zero polynomial."),
             ([0, 1], "P(x) may have zeros very close to stability boundary."),
             ([1, 0], "The degree of P(x) has been reduced to 0")]
        for P, m in T:
            with warnings.catch_warnings(record=True) as w:
                (dp, stable, nz) = mc01td('C', len(P)-1, P)
                self.assertEqual(str(w[0].message), m)


if __name__ == "__main__":
    unittest.main()
