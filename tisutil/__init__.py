from continuous import *
from tisutil.release import __version__

def use_python_types():
    mfu.nw.use_python_types()

def use_mpmath_types(dps=mfu.nw.dps_default_mpmath):
    mfu.nw.use_mpmath_types(dps)

def set_type_mode(mode, dps=None):
    mfu.nw.set_type_mode(mode, dps)
