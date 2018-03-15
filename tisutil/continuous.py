from discrete import *

class cSmat(mfu.cMat):
    def _getDiscreteContainer(self):
        return dSmat(units=self.units)

class cKmat(mfu.cMat):
    def _getDiscreteContainer(self):
        return dKmat(units=self.units)

class cTmat(mfu.cMat):
    def _getDiscreteContainer(self):
        return dTmat(units=self.units)

class cPolySmat(mfu.cPolyMat):
    def __init__(self, symMat, symVar, chanCalc):
        mfu.cMat.__init__(self,
          lambda ene: mfu.nw.fromSympyMatrix(symMat.subs(symVar, 
                                                         chanCalc.fk(ene))),
          chanCalc.getUnits())
        self.symMat = symMat
        self.symVar = symVar

    def _getDiscreteContainer(self):
        return dSmat(units=self.units)


def usePythonTypes_c(dps=mfu.nw.dps_default_python):
    mfu.usePythonTypes(dps)

def useMpmathTypes_c(dps=mfu.nw.dps_default_mpmath):
    mfu.useMpmathTypes(dps)
