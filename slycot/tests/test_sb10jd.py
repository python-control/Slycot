# ===================================================
# sb10jd tests

import unittest
from slycot import synthesis
import numpy as np
from numpy.testing import assert_raises, assert_almost_equal, assert_equal

# test1 input parameters

test1_n = 6
test1_m = 1
test1_np = 6

test1_A = np.array([[ 0,  0,  0, -1,  1,  0],
                    [ 0, 32,  0,  0, -1,  1],
                    [ 0,  0,  1,  0,  0,  0],
                    [ 0,  0,  0,  1,  0,  0],
                    [-1,  1,  0,  0,  0,  0],
                    [ 0, -1,  1,  0,  0,  0]])

    
test1_E = np.array([[  0,   0,   0,   0,   0,   0],
                    [  0,   0,   0,   0,   0,   0],
                    [  0,   0,   0, -10,   0,  10],
                    [  0,   0,   0,   0,   0,   0],
                    [  0,   0,   0,   0,   0,   0],
                    [  0,   0,   0,   0,   0,   0]])
                
test1_B = np.array([[-7.1],
                    [ 0. ],
                    [ 0. ],
                    [ 0. ],
                    [ 0. ],
                    [ 0. ]])
                
test1_C = np.eye(6)

test1_D = np.zeros((7,1))

# test1 expected results

test1_Aexp = np.array([[-0.00312500]])
test1_Bexp = np.array([[ 0.05899985]])
test1_Cexp = np.array([[-1.17518847e-02],   
                       [-1.17518847e-02],
                       [-1.17518847e-02],
                       [ 0.00000000e+00],
                       [ 0.00000000e+00],
                       [ 3.76060309e-01]])
test1_Dexp = np.array([[ 2.21875000e-01],
                       [ 2.21875000e-01],
                       [ 2.21875000e-01],
                       [ 0.00000000e+00],
                       [ 7.10000000e+00],
                       [ 0.00000000e+00]])

class test_sb10jd(unittest.TestCase):
    def test1_sb10jd(self):
        """ verify the output of sb10jd for a descriptor system """
        A,B,C,D = synthesis.sb10jd(test1_n,test1_m,test1_np,test1_A,test1_B,test1_C,test1_D,test1_E)
        assert_almost_equal(A, test1_Aexp)
        assert_almost_equal(B, test1_Bexp)
        assert_almost_equal(C, test1_Cexp)
        assert_almost_equal(D, test1_Dexp)

def suite():
   return unittest.TestLoader().loadTestsFromTestCase(TestConvert)


if __name__ == "__main__":
    unittest.main()
