import os
import sys
basedir = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,basedir+'/../..')

import dum_mats

a = dum_mats.row_offset_col_gain_zeroimag()
a.plot()

b = a.create_reduced_dim(0)
b.set_chart_title("Test title")
b.plot(logx=True, logy=True)

c = b.create_reduced_dim(1)
c.set_chart_parameters(leg_prefix="Test", xsize=10, ysize=10)
c.plot(imag=True)

d = a.trace()
d.plot()

e = a[1:].to_dKmat().eigenvalues()
e.plot()

f = a[1:].to_dXSmat().to_dXSsca()
f.plot()

g = a[1:].to_dXSsca()
g.plot()

h = a[1:].to_dEPhaseMat().to_dEPhaseSca()
h.plot()

i = a[1:].to_dEPhaseSca()
i.plot()

j = a[1:].to_dQmat()
j.plot()

k = j.eigenvalues()
k.plot()
