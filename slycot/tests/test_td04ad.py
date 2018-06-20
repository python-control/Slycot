#!/usr/bin/env python
#
# test_td04ad.py - test suite for tf -> ss conversion
# RvP, 04 Jun 2018

from __future__ import print_function

import unittest
from slycot import transform
import numpy as np

from numpy.testing import assert_raises, assert_almost_equal

class TestTf2SS(unittest.TestCase):

    def test_td04ad_case1(self):
        """td04ad: Convert with both 'C' and 'R' options"""
        
        # for octave:
        """
        num = { [0.0,  0.0, 1.0 ], [ 1.0, 0.0 ];
                [3.0, -1.0, 1.0 ], [ 0.0, 1.0 ];
                [0.0, 0.0, 1.0],   [ 0.0, 2.0 ] };
        den = { [1.0,  0.4, 3.0],  [ 1.0, 1.0 ];
                [1.0,  0.4, 3.0],  [ 1.0, 1.0 ];
                [1.0,  0.4, 3.0],  [ 1.0, 1.0 ]};
        """
        
        # common denominators for the inputs
        n = 2
        m = 2
        p = 3
        num = np.array([
            [ [0.0,  0.0, 1.0 ], [ 1.0, 0.0, 0.0 ] ],
            [ [3.0, -1.0, 1.0 ], [ 0.0, 1.0, 0.0 ] ],
            [ [0.0, 0.0, 1.0],   [ 0.0, 2.0, 0.0 ] ] ])
        p, m, d = num.shape
        numc = np.zeros((max(1, m, p), max(1, m, p), d), dtype=float)
        numc[:p,:m,:] = num

        denc = np.array(
            [ [1.0,  0.4, 3.0],  [ 1.0, 1.0, 0.0 ] ])
        indc = np.array(
            [ 2, 1 ], dtype=int)
        denr = np.array(
            [ [1.0,  0.4, 3.0],  [ 1.0, 1.0, 0.0 ], [1.0, 0.0, 0.0] ])
        indr = np.array(
            [ 2, 1, 0 ], dtype=int)

        n, A, B, C, D = transform.td04ad('C', 2, 3, indc, denc, numc)
        #print('A=\n', A, '\nB=\n', B, '\nC=\n', C, '\nD=\n', D)
        Ac = [ [-1, 0, 0], [ 0, -0.4, -0.3], [ 0, 10, 0]]
        Bc = [ [0, -1] ,[ 1 , 0], [ 0, 0]]
        Cc = [ [1, 0, 0.1], [-1, -2.2, -0.8], [ -2, 0, 0.1] ]
        Dc = [ [0, 1], [ 3, 0], [ 0, 0]]
        np.testing.assert_array_almost_equal(A, Ac)
        np.testing.assert_array_almost_equal(B, Bc)
        np.testing.assert_array_almost_equal(C, Cc)
        np.testing.assert_array_almost_equal(D, Dc)

        resr = transform.td04ad('R', 2, 3, indr, denr, num)
        #print('A=\n', A, '\nB=\n', B, '\nC=\n', C, '\nD=\n', D)

    def test_staticgain(self):
        """td04ad: Convert a transferfunction to SS with only static gain"""
        
        # 2 inputs, 3 outputs? columns share a denominator
        num = np.array([ [ [1.0], [2.0] ],
                         [ [0.2], [4.3] ],
                         [ [1.2], [3.2] ] ])
        p, m, d = num.shape
        numc = np.zeros((max(1, m, p), max(1, m, p), d), dtype=float)
        numc[:p,:m,:] = num
        
        # denc, columns share a common denominator
        denc = np.array([ [ 1.0], [0.5] ])
        Dc = (num / denc).reshape((3,2))
        idxc = np.zeros((2,), dtype=int)
        
        # denr, rows share a common denominator
        denr = np.array([ [1.0], [0.5], [3.0] ])
        idxr = np.zeros((3,), dtype=int)
        Dr = (num / denr[:, np.newaxis]).reshape((3,2))

        # fails with:
        # On entry to TB01XD parameter number  5 had an illegal value
        
        n, A, B, C, D = transform.td04ad('C', 2, 3, idxc, denc, numc)
        #print('A=\n', A, '\nB=\n', B, '\nC=\n', C, '\nD=\n', D)
        self.assertEqual(A.shape, (0,0))
        self.assertEqual(B.shape, (0,2))
        self.assertEqual(C.shape, (3,0))
        np.testing.assert_array_almost_equal(D, Dc)
        
        n, A, B, C, D = transform.td04ad('R', 2, 3, idxr, denr, num)
        #print('A=\n', A, '\nB=\n', B, '\nC=\n', C, '\nD=\n', D)
        self.assertEqual(A.shape, (0,0))
        self.assertEqual(B.shape, (0,2))
        self.assertEqual(C.shape, (3,0))
        np.testing.assert_array_almost_equal(D, Dr)
        

    def test_mixfeedthrough(self):
        """Test case popping up from control testing"""
        # a mix of feedthrough and dynamics. The problem from the control
        # package was somewhere else
        num = np.array([ [ [ 0.0,  0.0 ], [ 0.0, -0.2 ] ],
                         [ [ -0.1,  0.0 ], [ 0.0,  0.0 ] ] ])
        p, m, d = num.shape
        numc = np.zeros((max(1, m, p), max(1, m, p), d), dtype=float)
        numc[:p,:m,:] = num
        denc = np.array([ [ 1.0,  1.1 ], [ 1.0, 0.0 ] ])
        idxc = np.array([ 1, 0 ])
        n, A, B, C, D = transform.td04ad('C', 2, 2, idxc, denc, numc)
        np.testing.assert_array_almost_equal(D, np.array([[0,  0],[-0.1, 0]]))
        
    def test_toandfrom(self):

        A = np.array([[-3.0]])
        B = np.array([[0.1, 0.0]])
        C = np.array([[1.0],[0.0]])
        D = np.array([[0.0, 0.0],[0.0, 1.0]])

        tfout = transform.tb04ad(1, 2, 2, A, B, C, D)

        num = tfout[6]
        den = tfout[5]
        idxc = np.array([1, 0])
        n, At, Bt, Ct, Dt = transform.td04ad('R', 2, 2, idxc, den, num)
        np.testing.assert_array_almost_equal(D, Dt)
        np.testing.assert_array_almost_equal(A, At)
        
def suite():
   return unittest.TestLoader().loadTestsFromTestCase(TestTF2SS)


if __name__ == "__main__":
    unittest.main()
