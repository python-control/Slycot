# ab09nd - model order reduction

import numpy as np
from slycot import ab09nd

# SLICOT reference test; see SLICOT-Reference/examples/AB09ND.dat, AB09ND.res, TAB09ND.f
def test_slicot_ref():
    n = 7
    m = 2
    p = 3
    nr = None # Slycot uses None for ordsel = 'A'
    alpha = -0.6
    tol1 = 1e-1
    tol2 = 1e-14
    dico = 'C'
    job = 'N'
    equil = 'N'

    a = np.array([[-0.04165, 0.0000, 4.9200, -4.9200, 0.0000, 0.0000, 0.0000],
                  [-5.2100, -12.500, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000],
                  [0.0000, 3.3300, -3.3300, 0.0000, 0.0000, 0.0000, 0.0000],
                  [0.5450, 0.0000, 0.0000, 0.0000, -0.5450, 0.0000, 0.0000],
                  [0.0000, 0.0000, 0.0000, 4.9200, -0.04165, 0.0000, 4.9200],
                  [0.0000, 0.0000, 0.0000, 0.0000, -5.2100, -12.500, 0.0000],
                  [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 3.3300, -3.3300]])

    b = np.array([[0.0000, 0.0000],
                  [12.500, 0.0000],
                  [0.0000, 0.0000],
                  [0.0000, 0.0000],
                  [0.0000, 0.0000],
                  [0.0000, 12.500],
                  [0.0000, 0.0000]])

    c = np.array([[1.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000],
                  [0.0000, 0.0000, 0.0000, 1.0000, 0.0000, 0.0000, 0.0000],
                  [0.0000, 0.0000, 0.0000, 0.0000, 1.0000, 0.0000, 0.0000]])

    d = np.zeros((3,2))

    nr, ar, br, cr, dr, ns, hsv = \
            ab09nd(dico, job, equil, n, m, p, a, b, c, d, alpha, nr, tol1, tol2)

    # reference values
    ref_nr = 5
    ref_hsv = np.array([1.9178, 0.8621, 0.7666, 0.0336, 0.0246])
    ref_ar = np.array([[-0.5181, -1.1084,  0.0000,  0.0000,  0.0000],
                       [ 8.8157, -0.5181,  0.0000,  0.0000,  0.0000],
                       [ 0.0000,  0.0000,  0.5847,  0.0000,  1.9230],
                       [ 0.0000,  0.0000,  0.0000, -1.6606,  0.0000],
                       [ 0.0000,  0.0000, -4.3823,  0.0000, -3.2922]])

    ref_br = np.array([[-1.2837,  1.2837],
                       [-0.7522,  0.7522],
                       [-0.6379, -0.6379],
                       [ 2.0656, -2.0656],
                       [-3.9315, -3.9315]])

    ref_cr = np.array([[-0.1380, -0.6445, -0.6416, -0.6293, 0.2526],
                       [ 0.6246,  0.0196,  0.0000,  0.4107, 0.0000],
                       [ 0.1380,  0.6445, -0.6416,  0.6293, 0.2526]])

    ref_dr = np.array([[ 0.0582, -0.0090],
                       [ 0.0015, -0.0015],
                       [-0.0090,  0.0582]])

    assert nr == ref_nr

    np.testing.assert_array_almost_equal(hsv[:nr], ref_hsv, decimal=4)
    np.testing.assert_array_almost_equal(ar, ref_ar, decimal=4)
    np.testing.assert_array_almost_equal(br, ref_br, decimal=4)
    np.testing.assert_array_almost_equal(cr, ref_cr, decimal=4)
    np.testing.assert_array_almost_equal(dr, ref_dr, decimal=4)


# gh-242 regression test
# iwork was incorrectly sized
def test_gh242_regression():
    n = 67
    m = 1
    p = 1

    a = -np.eye(n)
    b = np.zeros((n, m))
    c = np.zeros((p, n))
    d = np.array([[42.24]])

    nr, ar, br, cr, dr, ns, hsv = \
        ab09nd(dico='C', job='B', equil='S', n=a.shape[0],
               m=b.shape[1], p=c.shape[0], A=a, B=b, C=c, D=d)

    assert nr == 0
    np.testing.assert_equal(d, dr)
