import os
import sys
basedir = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,basedir+'/../..')

import dum_mats

import unittest

nw = dum_mats.tu.mfu.nw

def use_mpmath_types():
    nw.use_mpmath_types()

class test_StoT(unittest.TestCase):
    def runTest(self):
        use_mpmath_types()
        d1 = dum_mats.row_offset_col_gain_zeroimag_Smat(sz=10)
        testdps = 1e-90
        for i in range(len(d1)):
            self.assertTrue(nw.are_matrices_close(d1[i][1],
                                                  d1.to_dTmat().to_dSmat()[i][1],
                                                  rtol=testdps, atol=testdps))

class test_StoK(unittest.TestCase):
    def runTest(self):
        use_mpmath_types()
        d1 = dum_mats.row_offset_col_gain_zeroimag_Smat(sz=10)
        testdps = 1e-90
        for i in range(len(d1)):
            self.assertTrue(nw.are_matrices_close(d1[i][1],
                                                  d1.to_dKmat().to_dSmat()[i][1],
                                                  rtol=testdps, atol=testdps))

class test_KtoS(unittest.TestCase):
    def runTest(self):
        use_mpmath_types()
        d1 = dum_mats.row_offset_col_gain_zeroimag_Kmat(sz=10)
        testdps = 1e-90
        for i in range(len(d1)):
            self.assertTrue(nw.are_matrices_close(d1[i][1],
                                                  d1.to_dSmat().to_dKmat()[i][1],
                                                  rtol=testdps, atol=testdps))

class test_KtoT(unittest.TestCase):
    def runTest(self):
        use_mpmath_types()
        d1 = dum_mats.row_offset_col_gain_zeroimag_Kmat(sz=10)
        testdps = 1e-90
        for i in range(len(d1)):
            self.assertTrue(nw.are_matrices_close(d1[i][1],
                                                  d1.to_dTmat().to_dKmat()[i][1],
                                                  rtol=testdps, atol=testdps))

class test_TtoS(unittest.TestCase):
    def runTest(self):
        use_mpmath_types()
        d1 = dum_mats.row_offset_col_gain_zeroimag_Tmat(sz=10)
        testdps = 1e-90
        for i in range(len(d1)):
            self.assertTrue(nw.are_matrices_close(d1[i][1],
                                                  d1.to_dSmat().to_dTmat()[i][1],
                                                  rtol=testdps, atol=testdps))

class test_TtoK(unittest.TestCase):
    def runTest(self):
        use_mpmath_types()
        d1 = dum_mats.row_offset_col_gain_posImag_Tmat(sz=5)
        testdps = 1e-90
        a = d1.to_dKmat()
        b = a.to_dTmat()
        for i in range(len(d1)):
            self.assertTrue(nw.are_matrices_close(d1[i][1],
                                                  d1.to_dKmat().to_dTmat()[i][1],
                                                  rtol=testdps, atol=testdps))

if __name__ == "__main__":
    #Just for debug
    b = test_TtoK()
    b.runTest()
