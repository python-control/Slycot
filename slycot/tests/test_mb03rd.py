#!/usr/bin/env python
#
# test_mb03rd.py - test suite for Shur form reduction
# RvP, 31 Jul 2019
import unittest
from slycot import transform
import numpy as np
from numpy.testing import assert_raises, assert_almost_equal, assert_equal
from scipy.linalg import schur

test1_A = np.array([
   [ 1.,  -1.,   1.,   2.,   3.,   1.,   2.,   3.],
   [ 1.,   1.,   3.,   4.,   2.,   3.,   4.,   2.],
   [ 0.,   0.,   1.,  -1.,   1.,   5.,   4.,   1.],
   [ 0.,   0.,   0.,   1.,  -1.,   3.,   1.,   2.],
   [ 0.,   0.,   0.,   1.,   1.,   2.,   3.,  -1.],
   [ 0.,   0.,   0.,   0.,   0.,   1.,   5.,   1.],
   [ 0.,   0.,   0.,   0.,   0.,   0.,   0.99999999, -0.99999999 ],
   [ 0.,   0.,   0.,   0.,   0.,   0.,   0.99999999,  0.99999999 ]
   ])
test1_n = test1_A.shape[0]

test1_Ar = np.array([
    [ 1.0000, -1.0000, -1.2247, -0.7071, -3.4186,  1.4577,  0.0000,  0.0000 ],
    [ 1.0000,  1.0000,  0.0000,  1.4142, -5.1390,  3.1637,  0.0000,  0.0000 ],
    [ 0.0000,  0.0000,  1.0000, -1.7321, -0.0016,  2.0701,  0.0000,  0.0000 ],
    [ 0.0000,  0.0000,  0.5774,  1.0000,  0.7516,  1.1379,  0.0000,  0.0000 ],
    [ 0.0000,  0.0000,  0.0000,  0.0000,  1.0000, -5.8606,  0.0000,  0.0000 ],
    [ 0.0000,  0.0000,  0.0000,  0.0000,  0.1706,  1.0000,  0.0000,  0.0000 ],
    [ 0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  1.0000, -0.8850 ],
    [ 0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  1.0000 ],
    ])

test1_Xr = np.array([
    [ 1.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.9045,  0.1957 ],
    [ 0.0000,  1.0000,  0.0000,  0.0000,  0.0000,  0.0000, -0.3015,  0.9755 ],
    [ 0.0000,  0.0000,  0.8165,  0.0000, -0.5768, -0.0156, -0.3015,  0.0148 ],
    [ 0.0000,  0.0000, -0.4082,  0.7071, -0.5768, -0.0156,  0.0000, -0.0534 ],
    [ 0.0000,  0.0000, -0.4082, -0.7071, -0.5768, -0.0156,  0.0000,  0.0801 ],
    [ 0.0000,  0.0000,  0.0000,  0.0000, -0.0276,  0.9805,  0.0000,  0.0267 ],
    [ 0.0000,  0.0000,  0.0000,  0.0000,  0.0332, -0.0066,  0.0000,  0.0000 ],
    [ 0.0000,  0.0000,  0.0000,  0.0000,  0.0011,  0.1948,  0.0000,  0.0000 ]
    ])

test1_pmax = 1e3
test1_tol = 0.01
class test_mb03rd(unittest.TestCase):
    def test1(self):
        # create schur form with scipy
        A, X = schur(test1_A)
        Ah, Xh = np.copy(A), np.copy(X)
        # on this basis, get the transform
        Ar, Xr, blks, eig = transform.mb03rd(
            test1_n, A, X, 'U', 'S', test1_pmax, test1_tol)
        # ensure X and A are unchanged
        assert_equal(A, Ah)
        assert_equal(X, Xh)
        # compare to test case results
        assert_almost_equal(Ar, test1_Ar, decimal=4)
        assert_almost_equal(Xr, test1_Xr, decimal=4)
        
if __name__ == "__main__":
    unittest.main()
