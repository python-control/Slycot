# ===================================================
# ag08bd tests

import unittest
from slycot import analysis
import numpy as np

from numpy.testing import assert_raises, assert_almost_equal, assert_equal

# test input parameters

test_A = np.array([[1, 0, 0, 0, 0, 0],
                   [0, 1, 0, 0, 0, 0],
                   [0, 0, 3, 0, 0, 0],
                   [0, 0, 0,-4, 0, 0],
                   [0, 0, 0, 0,-1, 0],
                   [0, 0, 0, 0, 1, 3]])

test_B = np.array([[0 , -1], 
                   [-1,  0],
                   [ 1, -1],
                   [ 0,  0],
                   [ 0,  1],
                   [-1, -1]])

test_C = np.array([[1, 0, 0, 1, 0, 0],
                   [0, 1, 0, 1, 0, 1],
                   [0, 0, 1, 0, 0, 1]])

test_D = np.zeros((3, 2))

test_A = test_A.astype(np.complex128)
test_B = test_B.astype(np.complex128)
test_C = test_C.astype(np.complex128)
test_D = test_D.astype(np.complex128)

                
class test_ab08n(unittest.TestCase):
    """ test1 to 4: Verify ag08bd with input parameters according to example in documentation """

    def test_ab08nd(self):
        #test [A-lambda*E]
        #B,C,D must have correct dimensions according to l,n,m and p, but cannot have zero length in any dimenstion. Then the wrapper will complain. The length is then set to one. 

        nu,rank,dinfz,nkror,nkrol,infz,kronr,kronl,Af,Bf = analysis.ab08nd(6,2,3,test_A,test_B,test_C,test_D)

    def test_ab08nz(self):
        #test [A-lambda*E]
        #B,C,D must have correct dimensions according to l,n,m and p, but cannot have zero length in any dimenstion. Then the wrapper will complain. The length is then set to one. 
        nu,rank,dinfz,nkror,nkrol,infz,kronr,kronl,Af,Bf = analysis.ab08nz(6,2,3,test_A,test_B,test_C,test_D)


def suite():
   return unittest.TestLoader().loadTestsFromTestCase(TestConvert)


if __name__ == "__main__":
    unittest.main()
