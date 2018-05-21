import os
import sys
basedir = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,basedir+'/../..')

import dum_mats

import unittest

nw = dum_mats.tu.mfu.nw

def use_mpmath_types():
    nw.use_mpmath_types()

class test_unitary(unittest.TestCase):
    def runTest(self):
        use_mpmath_types()
        d1 = dum_mats.unitary_matrix()
        testdps = 1e-90
        self.assertTrue(d1.to_dUniOpMat().is_unitary(rtol=testdps,
                                                     atol=testdps))
        self.assertTrue(nw.are_matrices_close(nw.identity(3),
                                              d1.to_dUniOpMat()[0][1],
                                              rtol=testdps, atol=testdps))

if __name__ == "__main__":
    #Just for debug
    b = test_unitary()
    b.runTest()
