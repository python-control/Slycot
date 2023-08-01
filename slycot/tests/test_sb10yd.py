import unittest
from slycot import synthesis
import numpy as np
from scipy import signal

from numpy.testing import assert_almost_equal, assert_equal

class test_sb10yd(unittest.TestCase):

    def test_sb10yd_exec(self):
        """Test execution. 
        """

        A = np.array([[0.0, 1.0], [-0.5, -0.1]])
        B = np.array([[0.0], [1.0]])
        C = np.array([[1.0, 0.0]])
        D = np.zeros((1,1))

        sys_tf = signal.ss2tf(A,B,C,D)
        num, den = sys_tf

        omega, H = signal.freqs(num.squeeze(), den)

        real_H_resp = np.real(H)
        imag_H_resp = np.imag(H)

        n = 2
        n_id, *_ = synthesis.sb10yd(
            0, 0, len(omega), 
            real_H_resp, imag_H_resp, omega, n, tol=0)

        np.testing.assert_equal(n, n_id)

    def test_sb10yd_allclose(self):
        """Compare given and identified frequency response.
        """

        A = np.array([[0.0, 1.0], [-0.5, -0.1]])
        B = np.array([[0.0], [1.0]])
        C = np.array([[1.0, 0.0]])
        D = np.zeros((1,1))

        sys_tf = signal.ss2tf(A,B,C,D)
        num, den = sys_tf

        omega, H = signal.freqs(num.squeeze(), den)

        real_H_resp = np.real(H)
        imag_H_resp = np.imag(H)

        n = 2
        n_id, A_id, B_id, C_id, D_id = synthesis.sb10yd(
            0, 0, len(omega), 
            real_H_resp, imag_H_resp, omega, n, tol=0)
        
        sys_tf_id = signal.ss2tf(A_id,B_id,C_id,D_id)
        num_id, den_id = sys_tf_id
        w_id, H_id = signal.freqs(num_id.squeeze(), den_id, worN=omega)

        np.testing.assert_allclose(abs(H),abs(H_id),rtol=0.3,atol=0)

if __name__ == "__main__":
    unittest.main()