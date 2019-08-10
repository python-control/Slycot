#!/usr/bin/env python
#
# test_mb.py - test suite for linear algebra commands
# bnavigator <code@bnavigator.de>, Aug 2019

import unittest
import numpy as np

from slycot import mb05md

from numpy.testing import assert_allclose


class test_mb(unittest.TestCase):

    def test_mb05md(self):
        """ test_mb05md: verify Matrix exponential with slicot doc example
        data from http://slicot.org/objects/software/shared/doc/MB05MD.html
        """
        A = np.array([[ 0.5,  0.,   2.3, -2.6],
                      [ 0.,   0.5, -1.4, -0.7],
                      [ 2.3, -1.4,  0.5,  0.0],
                      [-2.6, -0.7,  0.0,  0.5]])
        delta = 1.0
        Ar_ref = np.array([[  26.8551,  -3.2824,  18.7409, -19.4430],
                           [  -3.2824,   4.3474,  -5.1848,   0.2700],
                           [  18.7409,  -5.1848,  15.6012, -11.7228],
                           [ -19.4430,   0.2700, -11.7228,  15.6012]])
        Vr_ref = np.array([[-0.7,  0.7,  0.1, -0.1],
                           [ 0.1, -0.1,  0.7, -0.7],
                           [ 0.5,  0.5,  0.5,  0.5],
                           [-0.5, -0.5,  0.5,  0.5]])
        Yr_ref = np.array([[  -0.0349,   0.0050,   0.0249,  -0.0249],
                           [  38.2187,  -5.4598,  27.2991, -27.2991],
                           [   0.0368,   0.2575,   0.1839,   0.1839],
                           [  -0.7389,  -5.1723,   3.6945,   3.6945]])
        VAL_ref = np.array([-3., 4., -1., 2.])
        (Ar, Vr, Yr, VAL) = mb05md(A, delta)

        assert_allclose(Ar, Ar_ref, atol=0.0001)

        # Order of eigenvalues is not guaranteed, so we check them one by one.
        for i, e in enumerate(VAL):
            erow = np.ones(VAL.shape)*e
            i_ref = np.isclose(erow, VAL_ref)
            self.assertTrue(any(i_ref),
                            msg="eigenvalue {} not expetced".format(e))
            # Eigenvectors can have different scaling.
            vr_ref = Vr_ref[:, i_ref]*Vr[0, i]/Vr_ref[0, i_ref][0]
            assert_allclose(Vr[:, (i,)], vr_ref, atol=0.0001)

        assert_allclose(np.dot(Vr, Yr), np.dot(Vr_ref, Yr_ref), atol=0.0001)


if __name__ == "__main__":
    unittest.main()
