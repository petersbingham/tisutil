import os
import sys
basedir = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,basedir+'/../..')

import dum_mats
import tisutil as tu

import unittest

class test_convert_units(unittest.TestCase):
    def runTest(self):
        d_ryd = dum_mats.row_offset_col_gain_posNegImag()
        d_eV = d_ryd.convert_units(tu.eVs)
        for val in zip(d_ryd.sorted_keys(), d_eV.sorted_keys()):
            self.assertAlmostEqual(val[0]*tu.rydbergs_to_eVs, val[1])

if __name__ == "__main__":
    #Just for debug
    b = test_convert_units()
    b.runTest()
