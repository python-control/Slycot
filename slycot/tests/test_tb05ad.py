# ===================================================
# tb05ad tests
import unittest
from slycot import transform
import numpy as np

from numpy.testing import assert_raises, assert_almost_equal


# set the random seed so we can get consistent results.
np.random.seed(40)
CASES = {}

# This is a known failure for tb05ad when running job 'AG'
CASES['fail1'] = {'A': np.array([[-0.5,  0.,  0.,  0. ],
                              [ 0., -1.,  0. ,  0. ],
                              [ 1.,  0., -0.5,  0. ],
                              [ 0.,  1.,  0., -1. ]]),
                   'B': np.array([[ 1.,  0.],
                             [ 0.,  1.],
                             [ 0.,  0.],
                             [ 0.,  0.]]),
                   'C': np.array([[ 0.,  1.,  1.,  0.],
                             [ 0.,  1.,  0.,  1.],
                             [ 0.,  1.,  1.,  1.]])}

n = 20
p = 10
m = 14

CASES['pass1'] = {'A': np.random.randn(n, n),
                  'B':  np.random.randn(n, m),
                  'C':  np.random.randn(p, n)}


class test_tb05ad(unittest.TestCase):

    def test_tb05ad_ng(self):
        """
        Test that tb05ad with job 'NG' computes the correct
        frequency response.
        """
        for key in CASES:
            sys = CASES[key]
            self.check_tb05ad_AG_NG(sys, 10*1j, 'NG')

    @unittest.expectedFailure
    def test_tb05ad_ag_failure(self):
        """ Test tb05ad and job 'AG' (i.e., balancing enabled) fails
        on certain A matrices.
        """
        self.check_tb05ad_AG_NG(CASES['fail1'], 10*1j, 'AG')

    def test_tb05ad_nh(self):
        """Test that tb05ad with job = 'NH' computes the correct
        frequency response after conversion to Hessenberg form.

        First call tb05ad with job='NH' to transform to upper Hessenberg
        form which outputs the transformed system.
        Subsequently, call tb05ad with job='NH' using this transformed system.
        """
        jomega = 10*1j
        for key in CASES:
            sys = CASES[key]
            sys_transformed = self.check_tb05ad_AG_NG(sys, jomega, 'NG')
            self.check_tb05ad_NH(sys_transformed, sys, jomega)


    def test_tb05ad_errors(self):
        """
        Test tb05ad error handling. We give wrong inputs and
        and check that this raises an error.
        """
        self.check_tb05ad_errors(CASES['pass1'])

    def check_tb05ad_AG_NG(self, sys, jomega, job):
        """
        Check that tb05ad computes the correct frequency response when
        running jobs 'AG' and/or 'NG'.

        Inputs
        ------

        sys: A a dict of system matrices with keys 'A', 'B', and 'C'.
        jomega: A complex scalar, which is the frequency we are
                evaluating the system at.
        job: A string, either 'AG' or 'NH'

        Returns
        -------
        sys_transformed: A dict of the system matrices which have been
                         transformed according the job.
        """
        n, m = sys['B'].shape
        p = sys['C'].shape[0]
        result = transform.tb05ad(n, m, p, jomega,
                                  sys['A'], sys['B'], sys['C'], job=job)
        g_i = result[3]
        hinvb = np.linalg.solve(np.eye(n) * jomega - sys['A'], sys['B'])
        g_i_solve = sys['C'].dot(hinvb)
        assert_almost_equal(g_i_solve, g_i)
        sys_transformed = {'A': result[0], 'B': result[1], 'C': result[2]}
        return sys_transformed

    def check_tb05ad_NH(self, sys_transformed, sys, jomega):
        """
        Check tb05ad, computes the correct frequency response when
        job='NH' and we supply system matrices 'A', 'B', and 'C'
        which have been transformed by a previous call to tb05ad.
        We check we get the same result as computing C(sI - A)^-1B
        with the original system.

        Inputs
        ------

        sys_transformed: A a dict of the transformed (A in upper
                         hessenberg form) system matrices with keys
                         'A', 'B', and 'C'.

        sys: A dict of the original un-transformed system matrices.

        jomega: A complex scalar, which is the frequency to evaluate at.

        """

        n, m = sys_transformed['B'].shape
        p = sys_transformed['C'].shape[0]
        result = transform.tb05ad(n, m, p, jomega, sys_transformed['A'],
                                  sys_transformed['B'], sys_transformed['C'],
                                  job='NH')
        g_i = result[0]
        hinvb = np.linalg.solve(np.eye(n) * jomega - sys['A'], sys['B'])
        g_i_solve = sys['C'].dot(hinvb)
        assert_almost_equal(g_i_solve, g_i)

    def check_tb05ad_errors(self, sys):
        """
        Check the error handling of tb05ad. We give wrong inputs and
        and check that this raises an error.
        """
        n, m = sys['B'].shape
        p = sys['C'].shape[0]
        jomega = 10*1j
        # test error handling
        # wrong size A
        assert_raises(ValueError, transform.tb05ad, n+1, m, p,
                      jomega, sys['A'], sys['B'], sys['C'], job='NH')
        # wrong size B
        assert_raises(ValueError, transform.tb05ad, n, m+1, p,
                      jomega, sys['A'], sys['B'], sys['C'], job='NH')
        # wrong size C
        assert_raises(ValueError, transform.tb05ad, n, m, p+1,
                      jomega, sys['A'], sys['B'], sys['C'], job='NH')
        # unrecognized job
        assert_raises(ValueError, transform.tb05ad, n, m, p, jomega,
                      sys['A'], sys['B'], sys['C'], job='a')



def suite():
   return unittest.TestLoader().loadTestsFromTestCase(TestConvert)


if __name__ == "__main__":
    unittest.main()
