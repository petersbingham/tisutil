import os
import sys
from __builtin__ import raw_input
basedir = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,basedir+'/../..')

import dumMats

d = dumMats.row_offset_col_gain_zeroimag()
print str(d)

raw_input("Any key to continue.")
d = dumMats.row_offset_col_gain_posImag()
print str(d)

raw_input("Any key to continue.")
d = dumMats.row_offset_col_gain_negImag()
print str(d)

raw_input("Any key to continue.")
d = dumMats.row_offset_col_gain_posNegImag()
print str(d)