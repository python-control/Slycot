import numpy as np
from numpy.testing import assert_allclose

from slycot import analysis


class Test_ab04md:
    """Test ab04md.

    Example data taken from
    https://www.slicot.org/objects/software/shared/doc/AB04MD.html
    """

    Ac = np.array([[1.0, 0.5],
                   [0.5, 1.0]])
    Bc = np.array([[0.0, -1.0],
                   [1.0, 0.0]])
    Cc = np.array([[-1.0, 0.0],
                   [0.0, 1.0]])
    Dc = np.array([[1.0, 0.0],
                   [0.0, -1.0]])

    Ad = np.array([[-1.0, -4.0],
                   [-4.0, -1.0]])
    Bd = np.array([[2.8284, 0.0],
                   [0.0, -2.8284]])
    Cd = np.array([[0.0, 2.8284],
                   [-2.8284, 0.0]])
    Dd = np.array([[-1.0, 0.0],
                   [0.0, -3.0]])

    def test_ab04md_cont_disc_cont(self):
        """Test transformation from continuous - to discrete - to continuous time.
        """

        n, m = self.Bc.shape
        p = self.Cc.shape[0]

        Ad_t, Bd_t, Cd_t, Dd_t = analysis.ab04md(
            'C', n, m, p, self.Ac, self.Bc, self.Cc, self.Dc)

        Ac_t, Bc_t, Cc_t, Dc_t = analysis.ab04md(
            'D', n, m, p, Ad_t, Bd_t, Cd_t, Dd_t)

        assert_allclose(self.Ac, Ac_t)
        assert_allclose(self.Bc, Bc_t)
        assert_allclose(self.Cc, Cc_t)
        assert_allclose(self.Dc, Dc_t)

    def test_ab04md_disc_cont_disc(self):
        """Test transformation from discrete - to continuous - to discrete time.
        """

        n, m = self.Bc.shape
        p = self.Cc.shape[0]

        Ac_t, Bc_t, Cc_t, Dc_t = analysis.ab04md(
            'D', n, m, p, self.Ad, self.Bd, self.Cd, self.Dd)

        Ad_t, Bd_t, Cd_t, Dd_t = analysis.ab04md(
            'C', n, m, p, Ac_t, Bc_t, Cc_t, Dc_t)

        assert_allclose(self.Ad, Ad_t)
        assert_allclose(self.Bd, Bd_t)
        assert_allclose(self.Cd, Cd_t)
        assert_allclose(self.Dd, Dd_t)
