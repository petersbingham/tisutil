import os
import sys
basedir = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,basedir+'/../..')

import dumMats

import unittest

def useMpmathTypes():
    dumMats.tu.mfu.nw.useMpmathTypes()

def usePythonTypes():
    dumMats.tu.mfu.nw.usePythonTypes()

class test_StoT(unittest.TestCase):
    def runTest(self):
        useMpmathTypes()
        d1 = dumMats.rowOffsetColGain_zeroImag_Smat()
        d1.to_dTMats()
        usePythonTypes()

if __name__ == "__main__":
    #Just for debug
    b = test_StoT()
    b.runTest()
