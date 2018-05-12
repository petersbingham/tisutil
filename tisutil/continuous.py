from discrete import *

class cMat(mfu.cMat):
    def __init__(self, funPtr, asymcalc=None, sourceStr=""):
        mfu.cMat.__init__(self, funPtr, 
                          None if asymcalc is None else asymcalc.getUnits(), sourceStr)
        self.asymcalc = asymcalc

class cSmat(cMat):
    def _getDiscreteContainer(self):
        return dSmat(asymcalc=self.asymcalc)

class cKmat(cMat):
    def _getDiscreteContainer(self):
        return dKmat(asymcalc=self.asymcalc)

class cTmat(cMat):
    def _getDiscreteContainer(self):
        return dTmat(asymcalc=self.asymcalc)

def getContinuousScatteringMatrix(matType, funPtr, asymcalc, sourceStr=""):
    if matType == Smat:
        return cSmat(funPtr, asymcalc, sourceStr)
    elif matType == Kmat:
        return cKmat(funPtr, asymcalc, sourceStr)
    elif matType == Tmat:
        return cTmat(funPtr, asymcalc, sourceStr)
    else:
        raise Exception("Non-recognised matrix type.")

# k as in wavenumber, not K-matrix
class cPolykmat(mfu.cPolyMat):
    def __init__(self, symMat, symVar, asymcalc, sourceStr=""):
        mfu.cMat.__init__(self,
          lambda ene: mfu.nw.fromSympyMatrix(symMat.subs(symVar,
                                                         asymcalc.fk(ene))),
          asymcalc.getUnits(), sourceStr)
        self.asymcalc = asymcalc
        self.symMat = symMat
        self.symVar = symVar

class cPolySmat(cPolykmat):
    def _getDiscreteContainer(self):
        return dSmat(asymcalc=self.asymcalc)
