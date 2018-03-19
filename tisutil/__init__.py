from continuous import *

# These have to be repeated in the discrete and continuous containers since 
# both get their own copy of mfc.

def usePythonTypes(dps=mfu.nw.dps_default_python):
    usePythonTypes_d(dps)
    usePythonTypes_c(dps)

def useMpmathTypes(dps=mfu.nw.dps_default_mpmath):
    useMpmathTypes_d(dps)
    useMpmathTypes_c(dps)

def setTypeMode(mode, dps=None):
    setTypeMode_d(mode, dps)
    setTypeMode_c(mode, dps)
