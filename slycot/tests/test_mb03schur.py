#!/usr/bin/env python
"""
test_mb03schur.py
Created on Sun Jan 26 17:38:08 2020

@author: bnavigator

"""

import unittest
from slycot import math
import numpy as np

from numpy.testing import assert_allclose


class test_mb03schur(unittest.TestCase):
    """unit tests for schur decomposition functions"""

    def test_mb03vd_mb03vy_ex(self):
        """Test MB03VD and MB03VY
           with the example given in the MB03VD SLICOT documentation"""

        n = 4
        p = 2
        ilo = 1
        ihi = 4
        A = np.zeros((n, n, p))
        A[:, :, 0] = [[1.5, -.7, 3.5, -.7],
                      [1. , 0. , 2. , 3. ],
                      [1.5, -.7, 2.5, -.3],
                      [1. , 0. , 2. , 1. ]]
        A[:, :, 1] = [[1.5, -.7, 3.5, -.7],
                      [1. , 0. , 2. , 3. ],
                      [1.5, -.7, 2.5, -.3],
                      [1. , 0. , 2. , 1. ]]

        H_ref = np.zeros((n, n, p))
        H_ref[:, :, 0] = [[-2.3926,  2.7042, -0.9598, -1.2335],
                          [ 4.1417, -1.7046,  1.3001, -1.3120],
                          [ 0.0000, -1.6247, -0.2534,  1.6453],
                          [ 0.0000,  0.0000, -0.0169, -0.4451]]

        H_ref[:, :, 1] = [[-2.5495,  2.3402,  4.7021,  0.2329],
                          [ 0.0000,  1.9725, -0.2483, -2.3493],
                          [ 0.0000,  0.0000, -0.6290, -0.5975],
                          [ 0.0000,  0.0000,  0.0000, -0.4426]]

        Q_ref = np.zeros((n, n, p))
        Q_ref[:, :, 0] = [[ 1.0000,  0.0000,  0.0000,  0.0000],
                          [ 0.0000, -0.7103,  0.5504, -0.4388],
                          [ 0.0000, -0.4735, -0.8349, -0.2807],
                          [ 0.0000, -0.5209,  0.0084,  0.8536]]

        Q_ref[:, :, 1] = [[-0.5883,  0.2947,  0.7528, -0.0145],
                          [-0.3922, -0.8070,  0.0009, -0.4415],
                          [-0.5883,  0.4292, -0.6329, -0.2630],
                          [-0.3922, -0.2788, -0.1809,  0.8577]]

        HQ, Tau = math.mb03vd(n, ilo, ihi, A)

        H = np.zeros_like(HQ)
        Q = np.zeros_like(HQ)

        for k in range(p):
            Q[:, :, k] = np.tril(HQ[:, :, k])
            if k == 0:
                H[:, :, k] = np.triu(HQ[:n, :n, k], -1)
            elif k > 0:
                H[:, :, k] = np.triu(HQ[:n, :n, k])
            assert_allclose(H[:, :, k], H_ref[:, :, k], atol=1e-4)

        Qr = math.mb03vy(n, ilo, ihi, Q, Tau)

        for k in range(p):
            assert_allclose(Qr[:, :, k], Q_ref[:, :, k], atol=1e-4)

        # Computer Error: too machine dependent to test to reference value
        # SSQ_ref = 2.93760e-15
        # SSQ = 0.
        # for k in range(p):
        #     kp1 = k+1
        #     if kp1 > p-1:
        #         kp1 = 0
        #     P = Qr[:, :, k].T.dot(A[: ,: ,k]).dot(Qr[: ,: ,kp1]) - H[: ,: ,k]
        #     SSQ = np.sqrt(SSQ**2 + np.linalg.norm(P,'fro')**2)

    def test_mb03wd_ex(self):
        """Test MB03WD with the example given in the SLICOT documentation"""

        n = 4
        p = 2
        ilo = 1
        ihi = 4
        iloz = 1
        ihiz = 4
        job = 'S'
        compz = 'V'
        A = np.zeros((n, n, p))
        A[:, :, 0] = [[1.5, -.7, 3.5, -.7],
                      [1. , 0. , 2. , 3. ],
                      [1.5, -.7, 2.5, -.3],
                      [1. , 0. , 2. , 1. ]]
        A[:, :, 1] = [[1.5, -.7, 3.5, -.7],
                      [1. , 0. , 2. , 3. ],
                      [1.5, -.7, 2.5, -.3],
                      [1. , 0. , 2. , 1. ]]

        W_ref = np.array([6.449861+7.817717J,
                          6.449861-7.817717J,
                          0.091315+0.000000J,
                          0.208964+0.000000J])

        T_ref = np.zeros((n, n, p))
        T_ref[:, :, 0] = [[  2.2112,  4.3718, -2.3362,  0.8907],
                          [ -0.9179,  2.7688, -0.6570, -2.2426],
                          [  0.0000,  0.0000,  0.3022,  0.1932],
                          [  0.0000,  0.0000,  0.0000, -0.4571]]

        T_ref[:, :, 1] = [[  2.9169,  3.4539,  2.2016,  1.2367],
                          [  0.0000,  3.4745,  1.0209, -2.0720],
                          [  0.0000,  0.0000,  0.3022, -0.1932],
                          [  0.0000,  0.0000,  0.0000, -0.4571]]

        Z_ref = np.zeros((n, n, p))
        Z_ref[:, :, 0] = [[  0.3493,  0.6751, -0.6490,  0.0327],
                          [  0.7483, -0.4863, -0.1249, -0.4336],
                          [  0.2939,  0.5504,  0.7148, -0.3158],
                          [  0.4813, -0.0700,  0.2286,  0.8433]]


        Z_ref[:, :, 1] = [[  0.2372,  0.7221,  0.6490,  0.0327],
                          [  0.8163, -0.3608,  0.1249, -0.4336],
                          [  0.2025,  0.5902, -0.7148, -0.3158],
                          [  0.4863,  0.0076, -0.2286,  0.8433]]

        HQ, Tau = math.mb03vd(n, ilo, ihi, A)
        Q = math.mb03vy(n, ilo, ihi, HQ, Tau)
        T, Z, W = math.mb03wd(job, compz, n, ilo, ihi, iloz, ihiz, HQ, Q)

        # TODO (?)
        # isolate eigenvalues with math.mb03wx

        assert_allclose(W, W_ref, atol=1e-5)
        assert_allclose(T, T_ref, atol=1e-4)
        assert_allclose(Z, Z_ref, atol=1e-4)

        # Computer Error: too machine dependent to test to reference value
        # SSQ_ref = 7.18432D-15
        # SSQ = 0.
        # for k in range(p):
        #     kp1 = k+1
        #     if kp1 > p-1:
        #         kp1 = 0
        #     P = Zrr[:, :, k].T.dot(A[: ,: ,k]).dot(Zrr[: ,: ,kp1]) - Hrr[: ,: ,k]
        #     SSQ = np.sqrt(SSQ**2 + np.linalg.norm(P,'fro')**2)
