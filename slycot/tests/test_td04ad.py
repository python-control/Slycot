#!/usr/bin/env python
#
# test_td04ad.py - test suite for tf -> ss conversion
# RvP, 04 Jun 2018

from __future__ import print_function

import unittest
from slycot import transform
import numpy as np

from numpy.testing import assert_raises, assert_almost_equal

class test_td04ad(unittest.TestCase):

    def test_td04ad_case1(self):
        # common denominators for the inputs
        n = 2
        m = 2
        p = 3
        num = np.array([
            [ [0.0,  0.0, 1.0 ], [ 0.0, 1.0, 0.0 ] ],
            [ [3.0, -1.0, 1.0 ], [ 0.0, 0.0, 1.0 ] ],
            [ [0.0, 0.0, 1.0],   [ 0.0, 0.0, 2.0 ] ] ])
        numc = np.zeros((3, 3, 3),dtype=float)
        numc[:,:2,:] = num
        denc = np.array(
            [ [1.0,  0.4, 3.0],  [ 1.0, 1.0, 0.0 ] ])
        indc = np.array(
            [ 2, 1 ], dtype=int)
        denr = np.array(
            [ [1.0,  0.4, 3.0],  [ 1.0, 1.0, 0.0 ], [1.0, 0.0, 0.0] ])
        indr = np.array(
            [ 2, 1, 0 ], dtype=int)

        n, A, B, C, D = transform.td04ad('C', 2, 3, indc, denc, numc)
        print('A=\n', A, '\nB=\n', B, '\nC=\n', C, '\nD=\n', D)
        
        resr = transform.td04ad('R', 2, 3, indr, denr, num)
        print(resr)

def suite():
   return unittest.TestLoader().loadTestsFromTestCase(TestConvert)


if __name__ == "__main__":
    unittest.main()
