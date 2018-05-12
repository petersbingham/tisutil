import os
import sys
basedir = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,basedir+'/../..')

import dumMats

import unittest

def use_mpmath_types():
    dumMats.tu.mfu.nw.use_mpmath_types()

def use_python_types():
    dumMats.tu.mfu.nw.use_python_types()

class test_StoT(unittest.TestCase):
    def runTest(self):
        use_mpmath_types()
        d1 = dumMats.row_offset_col_gain_zeroimag_Smat()
        d1.to_dTmat()
        use_python_types()

if __name__ == "__main__":
    #Just for debug
    b = test_StoT()
    b.runTest()
