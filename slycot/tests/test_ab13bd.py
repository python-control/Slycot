# ===================================================
# ab08n* tests

import numpy as np
from numpy.testing import assert_allclose, assert_array_equal
from scipy import linalg, signal

from slycot import analysis


class Test_ab13bd:
    """ Test regular pencil construction ab08nX with input parameters
    according to example in documentation """

    A = np.array([[0.0, 1.0],[-0.5, -0.1]])
    B = np.array([[0.],[1.]])
    C = np.eye(2)
    D = np.zeros((2,1))

    Ad, Bd, Cd, Dd, dt = signal.cont2discrete((A, B, C, D), 0.1, method='zoh')

    def test_no_change_args_ccase(self):
        """ ab13md must not change its arguments. continuous system case.
        """

        acopy = self.A.copy()
        bcopy = self.B.copy()
        ccopy = self.C.copy()
        dcopy = self.D.copy()

        dico = 'C'
        jobn = 'H'

        n, m = self.B.shape
        p = self.C.shape[0]

        analysis.ab13bd(dico, jobn, n, m, p, self.A, self.B, self.C, self.D)
        assert_array_equal(self.A, acopy)
        assert_array_equal(self.B, bcopy)
        assert_array_equal(self.C, ccopy)
        assert_array_equal(self.D, dcopy)

    def test_no_change_args_dcase(self):
        """ ab13md must not change its arguments. discrete system case.
        """

        acopy = self.Ad.copy()
        bcopy = self.Bd.copy()
        ccopy = self.Cd.copy()
        dcopy = self.Dd.copy()

        dico = 'D'
        jobn = 'H'

        n, m = self.Bd.shape
        p = self.Cd.shape[0]

        analysis.ab13bd(dico, jobn, n, m, p, self.Ad, self.Bd, self.Cd, self.Dd)
        assert_array_equal(self.Ad, acopy)
        assert_array_equal(self.Bd, bcopy)
        assert_array_equal(self.Cd, ccopy)
        assert_array_equal(self.Dd, dcopy)

    def test_ab13bd_2norm_ccase(self):
        """ Compare ab13md with scipy solution (Lyapunov Equation). 
        continuous system case.
        """

        A = self.A
        B = self.B
        C = self.C
        D = self.D

        n, m = self.B.shape
        p = self.C.shape[0]

        dico = 'C'
        jobn = 'H'

        h2norm = analysis.ab13bd(dico, jobn, n, m, p, A, B, C, D)

        Lc = linalg.solve_continuous_lyapunov(A, -B@B.T)    
        h2norm_Lc = np.sqrt(np.trace(C@Lc@C.T))
        print(h2norm_Lc, h2norm)
        assert_allclose(h2norm_Lc, h2norm, atol=1e-5)

        Lo = linalg.solve_continuous_lyapunov(A.T, -C.T@C)    
        h2norm_Lo = np.sqrt(np.trace(B.T@Lo@B))
        print(h2norm_Lo, h2norm)
        assert_allclose(h2norm_Lo, h2norm, atol=1e-5)

    def test_ab13bd_2norm_dcase(self):
        """ Compare ab13md with scipy solution (Lyapunov Equation). 
        discrete system case.
        """

        Ad = self.Ad
        Bd = self.Bd
        Cd = self.Cd
        Dd = self.Dd

        n, m = Bd.shape
        p = Cd.shape[0]

        dico = 'D'
        jobn = 'H'

        h2norm = analysis.ab13bd(dico, jobn, n, m, p, Ad, Bd, Cd, Dd)

        Lc = linalg.solve_discrete_lyapunov(Ad, Bd@Bd.T)    
        h2norm_Lc = np.sqrt(np.trace(Cd@Lc@Cd.T + Dd@Dd.T))
        print(h2norm, h2norm_Lc)
        assert_allclose(h2norm_Lc, h2norm, atol=1e-5)

        Lo = linalg.solve_discrete_lyapunov(Ad.T, Cd.T@Cd)    
        h2norm_Lo = np.sqrt(np.trace(Bd.T@Lo@Bd + Dd.T@Dd))
        print(h2norm, h2norm_Lo)
        assert_allclose(h2norm_Lo, h2norm, atol=1e-5)

