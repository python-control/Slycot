#
# test_mb.py - test suite for linear algebra commands
# bnavigator <code@bnavigator.de>, Aug 2019

import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.linalg import schur

from slycot import math, mb02ed, mb03rd, mb03vd, mb03vy, mb03wd, mb05md, mb05nd
from slycot.exceptions import SlycotArithmeticError, SlycotResultWarning, SlycotParameterError

from .test_exceptions import assert_docstring_parse


def test_mb02ed_ex():
    """Test MB02ED using the example given in the MB02ED SLICOT Documentation"""
    n = 3
    k = 3
    nrhs = 2
    TYPET = "C"
    T = np.array(
        [
            [3.0000, 1.0000, 0.2000],
            [1.0000, 4.0000, 0.4000],
            [0.2000, 0.4000, 5.0000],
            [0.1000, 0.1000, 0.2000],
            [0.2000, 0.0400, 0.0300],
            [0.0500, 0.2000, 0.1000],
            [0.1000, 0.0300, 0.1000],
            [0.0400, 0.0200, 0.2000],
            [0.0100, 0.0300, 0.0200],
        ]
    )
    B = np.array(
        [
            [1.0000, 2.0000],
            [1.0000, 2.0000],
            [1.0000, 2.0000],
            [1.0000, 2.0000],
            [1.0000, 2.0000],
            [1.0000, 2.0000],
            [1.0000, 2.0000],
            [1.0000, 2.0000],
            [1.0000, 2.0000],
        ]
    )
    X = np.array(
        [
            [0.2408, 0.4816],
            [0.1558, 0.3116],
            [0.1534, 0.3068],
            [0.2302, 0.4603],
            [0.1467, 0.2934],
            [0.1537, 0.3075],
            [0.2349, 0.4698],
            [0.1498, 0.2995],
            [0.1653, 0.3307],
        ]
    )

    result,_ = mb02ed(T=T, B=B, n=n, k=k, typet=TYPET, nrhs=nrhs)
    np.testing.assert_almost_equal(result, X, decimal=4)

def test_mb02ed_parameter_errors():
    """Test for errors in the input parameters of MB02ED"""
    n = 3
    k = 3
    nrhs = 2
    TYPET = "C"
    T = np.array(
        [
            [3.0000, 1.0000, 0.2000],
            [1.0000, 4.0000, 0.4000],
            [0.2000, 0.4000, 5.0000],
            [0.1000, 0.1000, 0.2000],
            [0.2000, 0.0400, 0.0300],
            [0.0500, 0.2000, 0.1000],
            [0.1000, 0.0300, 0.1000],
            [0.0400, 0.0200, 0.2000],
            [0.0100, 0.0300, 0.0200],
        ]
    )
    B = np.array(
        [
            [1.0000, 2.0000],
            [1.0000, 2.0000],
            [1.0000, 2.0000],
            [1.0000, 2.0000],
            [1.0000, 2.0000],
            [1.0000, 2.0000],
            [1.0000, 2.0000],
            [1.0000, 2.0000],
            [1.0000, 2.0000],
        ]
    )

    # Test for wrong parameter typet
    with pytest.raises(expected_exception=SlycotParameterError, match='typet must be either "R" or "C"') as cm:
        mb02ed(T=T, B=B, n=n, k=k, typet='U', nrhs=nrhs)
        assert cm.value.info == -1
    #Test for negative number of columns
    with pytest.raises(expected_exception=SlycotParameterError, match="k must be >= 0") as cm:
        mb02ed(T=T, B=B, n=n, k=-1, typet=TYPET, nrhs=nrhs)
        assert cm.value.info == -2
    #Test for negative number of blocks
    with pytest.raises(expected_exception=SlycotParameterError, match="n must be >= 0") as cm:
        mb02ed(T=T, B=B, n=-1, k=k, typet=TYPET, nrhs=nrhs)
        assert cm.value.info == -3
    #Test for negative number of right hand sides
    with pytest.raises(expected_exception=SlycotParameterError, match="nrhs must be >= 0") as cm:
        mb02ed(T=T, B=B, n=n, k=k, typet=TYPET, nrhs=-1)
        assert cm.value.info == -4


def test_mb02ed_matrix_error():
    """Test for a negative definite input matrix in MB02ED"""
    n = 3
    k = 3
    nrhs = 2
    TYPET = "C"
    T = np.array(
        [
            [3.0000, 1.0000, 0.2000],
            [1.0000, 4.0000, 0.4000],
            [0.2000, 0.4000, 5.0000],
            [0.1000, 0.1000, 0.2000],
            [0.2000, 0.0400, 0.0300],
            [0.0500, 0.2000, 0.1000],
            [0.1000, 0.0300, 0.1000],
            [0.0400, 0.0200, 0.2000],
            [0.0100, 0.0300, 0.0200],
        ]
    )
    # Create a negative definite matrix
    T = -1 * T
    B = np.array(
        [
            [1.0000, 2.0000],
            [1.0000, 2.0000],
            [1.0000, 2.0000],
            [1.0000, 2.0000],
            [1.0000, 2.0000],
            [1.0000, 2.0000],
            [1.0000, 2.0000],
            [1.0000, 2.0000],
            [1.0000, 2.0000],
        ]
    )

    with pytest.raises(SlycotArithmeticError,
                       match = "The reduction algorithm failed. "
                       "The Toeplitz matrix associated\nwith T "
                       r"is not numerically positive definite.") as cm:
        mb02ed(T=T, B=B, n=n, k=k, typet=TYPET, nrhs=nrhs)
        assert cm.value.info == 1


def test_mb03rd():
    """ Test for Schur form reduction.

    RvP, 31 Jul 2019"""

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

    test1_W = np.array([1+1j, 1-1j,
                        1+1j, 1-1j,
                        0.99999+0.99999j, 0.99999-0.99999j,
                        1., 1.])

    test1_pmax = 1e3
    test1_tol = 0.01
    # create schur form with scipy
    A, X = schur(test1_A)
    Ah, Xh = np.copy(A), np.copy(X)
    # on this basis, get the transform
    Ar, Xr, blsize, W = mb03rd(
            test1_n, A, X, 'U', 'S', test1_pmax, test1_tol)
    # ensure X and A are unchanged
    assert_allclose(A, Ah)
    assert_allclose(X, Xh)
    # compare to test case results
    assert_allclose(Ar, test1_Ar, atol=0.0001)
    assert_allclose(Xr, test1_Xr, atol=0.0001)
    assert_allclose(W,  test1_W, atol=0.0001)

    # Test that the non sorting options do not throw errors and that Xr is
    # returned as None for jobx='N'
    for sort in ['N', 'C', 'B']:
        Ar, Xr, blsize, W = mb03rd(
                test1_n, A, X, 'N', sort, test1_pmax, test1_tol)
        assert Xr is None

def test_mb03rd_default():
    # regression: mb03rd was failing with no third arg (X) supplied
    A = np.array([[ 6, -1, -7, -2,  2],
                    [-3,  4,  2, -7,  6],
                    [-6, -9, -3, -1, 10],
                    [-2, -4,  1,  5,  7],
                    [-7, -5, -6,  6,  7]])

    Aschur, Tschur = schur(A)

    X = Tschur.copy()

    Ar, Xr, blsize, W = mb03rd(Aschur.shape[0], Aschur, X, 'U', 'N', pmax=1.0, tol=0.0)

    Ar2, Xr2, blsize2, W2 = mb03rd(Aschur.shape[0], Aschur)

    assert_allclose(Ar, Ar2)
    assert_allclose(Xr, Tschur.dot(Xr2))

def test_mb03vd_mb03vy_ex():
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

    HQ, Tau = mb03vd(n, ilo, ihi, A)

    H = np.zeros_like(HQ)
    Q = np.zeros_like(HQ)

    for k in range(p):
        Q[:, :, k] = np.tril(HQ[:, :, k])
        if k == 0:
            H[:, :, k] = np.triu(HQ[:n, :n, k], -1)
        elif k > 0:
            H[:, :, k] = np.triu(HQ[:n, :n, k])
        assert_allclose(H[:, :, k], H_ref[:, :, k], atol=1e-4)

    Qr = mb03vy(n, ilo, ihi, Q, Tau)

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

def test_mb03wd_ex():
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

    HQ, Tau = mb03vd(n, ilo, ihi, A)
    Q = mb03vy(n, ilo, ihi, HQ, Tau)
    T, Z, W = mb03wd(job, compz, n, ilo, ihi, iloz, ihiz, HQ, Q)

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


def test_mb05md():
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
        assert any(i_ref), f"eigenvalue {e} not expected"
        # Eigenvectors can have different scaling.
        vr_ref = Vr_ref[:, i_ref]*Vr[0, i]/Vr_ref[0, i_ref][0]
        assert_allclose(Vr[:, (i,)], vr_ref, atol=0.0001)

    assert_allclose(np.dot(Vr, Yr), np.dot(Vr_ref, Yr_ref), atol=0.0001)

def test_mb05md_warning():
    """Check that the correct warning is raised from docstring"""
    A = np.diag([3., 3., 3., 3.]) + np.diag([1., 1., 1.], k=1)
    delta = 0.1

    with pytest.warns(SlycotResultWarning,
                        match="\n"
                            "Matrix A is defective, possibly "
                            "due to rounding errors.") as record:
        (Ar, Vr, Yr, VAL) = mb05md(A, delta)
    assert record[0].message.info == 6

def test_mb05nd():
    """ test_mb05nd: verify Matrix exponential and integral
    data from http://slicot.org/objects/software/shared/doc/MB05ND.html
    """
    A = np.array([[5.0,   4.0,   3.0,   2.0,   1.0],
                    [1.0,   6.0,   0.0,   4.0,   3.0],
                    [2.0,   0.0,   7.0,   6.0,   5.0],
                    [1.0,   3.0,   1.0,   8.0,   7.0],
                    [2.0,   5.0,   7.0,   1.0,   9.0]])
    delta = 0.1
    F_ref = np.array([[1.8391,  0.9476,  0.7920,  0.8216,  0.7811],
                        [0.3359,  2.2262,  0.4013,  1.0078,  1.0957],
                        [0.6335,  0.6776,  2.6933,  1.6155,  1.8502],
                        [0.4804,  1.1561,  0.9110,  2.7461,  2.0854],
                        [0.7105,  1.4244,  1.8835,  1.0966,  3.4134]])
    H_ref = np.array([[0.1347,  0.0352,  0.0284,  0.0272,  0.0231],
                        [0.0114,  0.1477,  0.0104,  0.0369,  0.0368],
                        [0.0218,  0.0178,  0.1624,  0.0580,  0.0619],
                        [0.0152,  0.0385,  0.0267,  0.1660,  0.0732],
                        [0.0240,  0.0503,  0.0679,  0.0317,  0.1863]])

    (F, H) = mb05nd(A, delta)

    assert_allclose(F, F_ref, atol=0.0001)
    assert_allclose(H, H_ref, atol=0.0001)


@pytest.mark.parametrize(
    'fun,          exception_class,       erange,         checkvars',
    ((math.mb03wd, SlycotResultWarning,   (4,),           {'ilo': 2,
                                                           'ihi': 5}),
     (math.mb05md, SlycotResultWarning,   (2, 3, 4),      {'n': 2}),
     (math.mb05nd, SlycotArithmeticError, (2, 3),         {'n': 2,
                                                           'delta': 1000}),
     (math.mc01td, SlycotResultWarning,   (1, 2, (1, 0)), {'dico': 'C'})))
def test_mb_docparse(fun, exception_class, erange, checkvars):
    assert_docstring_parse(fun.__doc__,  exception_class, erange, checkvars)

