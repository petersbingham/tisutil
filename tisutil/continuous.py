from discrete import *

class cMat(mfu.cMat):
    def __init__(self, funPtr, asymCal=None):
        mfu.cMat.__init__(self, funPtr, 
                          None if asymCal is None else asymCal.getUnits())
        self.asymCal = asymCal

class cSmat(cMat):
    def _getDiscreteContainer(self):
        return dSmat(asymCal=self.asymCal)

class cKmat(cMat):
    def _getDiscreteContainer(self):
        return dKmat(asymCal=self.asymCal)

class cTmat(cMat):
    def _getDiscreteContainer(self):
        return dTmat(asymCal=self.asymCal)

def getContinuousScatteringMatrix(matType, funPtr, asymCal):
    if matType == Smat:
        return cSmat(funPtr, asymCal)
    elif matType == Kmat:
        return cKmat(funPtr, asymCal)
    elif matType == Tmat:
        return cTmat(funPtr, asymCal)
    else:
        raise Exception("Non-recognised matrix type.")

# k as in wavenumber, not K-matrix
class cPolykmat(mfu.cPolyMat):
    def __init__(self, symMat, symVar, asymCal):
        mfu.cMat.__init__(self,
          lambda ene: mfu.nw.fromSympyMatrix(symMat.subs(symVar,
                                                         asymCal.fk(ene))),
          asymCal.getUnits())
        self.asymCal = asymCal
        self.symMat = symMat
        self.symVar = symVar

class cPolySmat(cPolykmat):
    def _getDiscreteContainer(self):
        return dSmat(asymCal=self.asymCal)
