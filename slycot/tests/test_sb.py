# ===================================================
# sb* synthesis tests

from slycot import synthesis
from numpy import array, eye, zeros
from numpy.testing import assert_allclose


def test_sb02mt():
    """Test if sb02mt is callable

    This is a dummy test, not really checking the wrapper of the FORTRAN
    function
    """
    out = synthesis.sb02mt(1, 1, 1., 1.)
    assert(len(out) == 8)


def test_sb10ad():
    """Test sb10ad, Hinf synthesis"""
    a = array([[-1]])
    b = array([[1, 1]])
    c = array([[1],
               [1]])
    d = array([[0, 1],
               [1, 0]])

    n = 1
    m = 2
    np = 2
    ncon = 1
    nmeas = 1
    gamma = 10

    gamma_est, Ak, Bk, Ck, Dk, Ac, Bc, Cc, Dc, rcond = synthesis.sb10ad(
        n, m, np, ncon, nmeas, gamma, a, b, c, d)
    # from Octave, which also uses SB10AD:
    #   a= -1; b1= 1; b2= 1; c1= 1; c2= 1; d11= 0; d12= 1; d21= 1; d22= 0;
    #   g = ss(a,[b1,b2],[c1;c2],[d11,d12;d21,d22]);
    #   [k,cl] = hinfsyn(g,1,1);
    # k.a is Ak, cl.a is Ac
    # gamma values don't match; not sure that's critical
    # this is a bit fragile
    # a simpler, more robust check might be to check stability of Ac
    assert_allclose(Ak, array([[-3.]]))
    assert_allclose(Ac, array([[-1., -1.],
                               [1., -3.]]))


def test_sb10jd():
    """ verify the output of sb10jd for a descriptor system """

    # test1 input parameters
    n = 6
    m = 1
    np = 6

    A = array([[ 0,  0,  0, -1,  1,  0],
               [ 0, 32,  0,  0, -1,  1],
               [ 0,  0,  1,  0,  0,  0],
               [ 0,  0,  0,  1,  0,  0],
               [-1,  1,  0,  0,  0,  0],
               [ 0, -1,  1,  0,  0,  0]])
    E = array([[  0,   0,   0,   0,   0,   0],
               [  0,   0,   0,   0,   0,   0],
               [  0,   0,   0, -10,   0,  10],
               [  0,   0,   0,   0,   0,   0],
               [  0,   0,   0,   0,   0,   0],
               [  0,   0,   0,   0,   0,   0]])
    B = array([[-7.1],
               [ 0. ],
               [ 0. ],
               [ 0. ],
               [ 0. ],
               [ 0. ]])
    C = eye(6)
    D = zeros((7, 1))

    # test1 expected results
    Aexp = array([[-0.003125]])
    Bexp = array([[ 0.059000]])
    Cexp = array([[-1.17519e-02],
                  [-1.17519e-02],
                  [-1.17519e-02],
                  [ 0.         ],
                  [ 0.         ],
                  [ 3.76060e-01]])
    Dexp = array([[ 2.21875e-01],
                  [ 2.21875e-01],
                  [ 2.21875e-01],
                  [ 0.         ],
                  [ 7.100000+00],
                  [ 0.         ]])

    A_r, B_r, C_r, D_r = synthesis.sb10jd(n, m, np, A, B, C, D, E)
    assert_allclose(A_r, Aexp, atol=1e-5)
    assert_allclose(B_r, Bexp, atol=1e-5)
    assert_allclose(C_r, Cexp, atol=1e-5)
    assert_allclose(D_r, Dexp, atol=1e-5)
