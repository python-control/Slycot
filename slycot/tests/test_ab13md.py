import numpy as np
import pytest
from numpy.testing import assert_allclose

from slycot import ab13md

# References:
# [1] Skogestand & Postlethwaite, Multivariable Feedback Control, 1996
# [2] slycot/src/SLICOT-Reference/examples

def slicot_example():
    # [2] AB13MD.dat & TAB13MD.f
    n = 6
    nblock = np.array([1, 1, 2, 1, 1])
    itype = np.array([1, 1, 2, 2, 2])
    # this unpleasant looking array is the result of text editor
    # search-and-replace on the Z array in AB13MD.dat
    Z = np.array([
        complex(-1.0e0,6.0e0),  complex(2.0e0,-3.0e0),  complex(3.0e0,8.0e0),
        complex(3.0e0,8.0e0),   complex(-5.0e0,-9.0e0), complex(-6.0e0,2.0e0),
        complex(4.0e0,2.0e0),   complex(-2.0e0,5.0e0),  complex(-6.0e0,-7.0e0),
        complex(-4.0e0,11.0e0), complex(8.0e0,-7.0e0),  complex(12.0e0,-1.0e0),
        complex(5.0e0,-4.0e0),  complex(-4.0e0,-8.0e0), complex(1.0e0,-3.0e0),
        complex(-6.0e0,14.0e0), complex(2.0e0,-5.0e0),  complex(4.0e0,16.0e0),
        complex(-1.0e0,6.0e0),  complex(2.0e0,-3.0e0),  complex(3.0e0,8.0e0),
        complex(3.0e0,8.0e0),   complex(-5.0e0,-9.0e0), complex(-6.0e0,2.0e0),
        complex(4.0e0,2.0e0),   complex(-2.0e0,5.0e0),  complex(-6.0e0,-7.0e0),
        complex(-4.0e0,11.0e0), complex(8.0e0,-7.0e0),  complex(12.0e0,-1.0e0),
        complex(5.0e0,-4.0e0),  complex(-4.0e0,-8.0e0), complex(1.0e0,-3.0e0),
        complex(-6.0e0,14.0e0), complex(2.0e0,-5.0e0),  complex(4.0e0,16.0e0),
    ])
    Z = np.reshape(Z, (n, n))
    return Z, nblock, itype


def test_cached_inputoutput():
    # check x, cached working area, input and output, and error
    Z, nblock, itype = slicot_example()

    m = len(nblock)
    mr = np.count_nonzero(1==itype)

    mu0, d0, g0, x0 = ab13md(Z, nblock, itype)
    assert m+mr-1 == len(x0)

    mu1, d1, g1, x1 = ab13md(Z, nblock, itype, x0)

    assert_allclose(mu1, mu0)

    with pytest.raises(ValueError):
        mu0, d, g, x = ab13md(Z, nblock, itype, np.ones(mr+mr))


class TestReference:
    # check a few reference cases
    def test_complex_scalar(self):
        # [1] (8.74)
        nblock = np.array([1])
        itype = np.array([2]) # complex

        z = np.array([[1+2j]])
        mu = ab13md(z,nblock,itype)[0]
        assert_allclose(mu, abs(z))


    def test_real_scalar_real_uncertainty(self):
        # [1] (8.75)
        nblock=np.array([1])
        itype=np.array([1]) # real
        z = np.array([[5.34]])
        mu = ab13md(z,nblock,itype)[0]
        assert_allclose(mu, abs(z))


    def test_complex_scalar_real_uncertainty(self):
        # [1] (8.75)
        nblock=np.array([1])
        itype=np.array([1]) # real

        z = np.array([[6.78j]])
        mu = ab13md(z,nblock,itype)[0]
        assert_allclose(mu, 0)


    def test_sp85_part1(self):
        # [1] Example 8.5 part 1, unstructured uncertainty
        M = np.array([[2, 2], [-1, -1]])
        muref = 3.162

        nblock=np.array([2])
        itype=np.array([2])

        mu = ab13md(M, nblock, itype)[0]
        assert_allclose(mu, muref, rtol=5e-4)


    def test_sp85_part2(self):
        # [1] Example 8.5 part 2, structured uncertainty
        M = np.array([[2, 2], [-1, -1]])
        muref = 3.000

        nblock=np.array([1, 1])
        itype=np.array([2, 2])

        mu = ab13md(M, nblock, itype)[0]
        assert_allclose(mu, muref, rtol=5e-4)

    def test_slicot(self):
        # besides muref, check that we've output d and g correctly

        muref =  0.4174753408e+02 # [2] AB13MD.res

        Z, nblock, itype = slicot_example()

        mu, d, g, x = ab13md(Z, nblock, itype)

        assert_allclose(mu, muref)

        ZH = Z.T.conj()
        D = np.diag(d)
        G = np.diag(g)

        # this matrix should be negative semi-definite
        negsemidef = (ZH @ D**2 @ Z
                      + 1j * (G@Z - ZH@G)
                      - mu**2 * D**2)

        # a check on the algebra, not ab13md!
        assert_allclose(negsemidef, negsemidef.T.conj())

        evals = np.linalg.eigvalsh(negsemidef)
        assert max(evals) < np.finfo(float).eps**0.5
