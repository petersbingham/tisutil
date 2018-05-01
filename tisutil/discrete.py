import matfuncutil as mfu
from channelutil.units import *
import copy

class dBase:
    def convertUnits(self, newUnits):
        if newUnits == self.units:
            return self
        elif self.units==RYDs:
            if newUnits==HARTs:
                fac = RYD_to_HART
            elif newUnits==eVs:
                fac = RYD_to_EV
            else:
                raise Exception("Unknown conversion")
        elif self.units==HARTs:
            if newUnits==RYDs:
                fac = 1./RYD_to_HART
            elif newUnits==eVs:
                fac = HART_to_EV
            else:
                raise Exception("Unknown conversion")
        elif self.units==eVs:
            if newUnits==RYDs:
                fac = 1./RYD_to_EV
            elif newUnits==HARTs:
                fac = 1./HART_to_EV
            else:
                raise Exception("Unknown conversion")
        else:
            raise Exception("Unknown conversion")
        newAsymCalc = copy.deepcopy(self.asymCal)
        newAsymCalc.units = newUnits
        newItem = self._createNewItem(newAsymCalc)
        self._initNewItem(newItem)
        for key,val in self.iteritems():
            newItem[key*fac] = val
        return newItem

class dVal(mfu.dVal, dBase):
    pass

class dVec(mfu.dVec, dBase):
    def _getReductionContainer(self):
        return dVal(units=self.units)

class dMat(mfu.dMat, dBase):
    def __init__(self, d={}, asymCal=None, sourceStr=""):
        mfu.dMat.__init__(self, d, 
                          None if asymCal is None else asymCal.getUnits(), sourceStr)
        self.asymCal = asymCal
    def getCheckStr(self):
        return mfu.dMat.getCheckStr(self) + "\n" + str(self.asymCal)
    def _getReductionContainer(self):
        return dVec(units=self.units)
    def _createNewItem(self, asymCal=None, newType=None):
        if asymCal is None:
            asymCal = self.asymCal
        if newType is None:
            newType = type(self)
        newItem = newType(asymCal=asymCal)
        newItem.sourceStr = self.sourceStr
        return newItem

class dSmat(dMat):
    def __init__(self, d={}, asymCal=None, sourceStr=""):
        dMat.__init__(self, d, asymCal, sourceStr)
        self.chartTitle = "S matrix"

    def to_dSmat(self):
        return self
    def to_dTmat(self):
        newItem = self._createNewItem(self.asymCal, newType=dTmat)
        self._initNewItem(newItem)
        for key in self:
            val = self[key] # force fun eval if relevant
            newItem[key] = val - mfu.nw.identity(mfu.nw.shape(val)[0])
        return newItem
    def to_dKmat(self):
        raise NotImplementedError

    def to_dXSmat(self):
        raise NotImplementedError
    def to_dEPhaseMat(self):
        raise NotImplementedError
    def to_dUniOpMat(self):
        raise NotImplementedError

class dKmat(dMat):
    def __init__(self, d={}, asymCal=None, sourceStr=""):
        dMat.__init__(self, d, asymCal, sourceStr)
        self.chartTitle = "K matrix"

    def to_dSmat(self):
        newItem = self._createNewItem(self.asymCal, newType=dSmat)
        self._initNewItem(newItem)
        for key in self:
            val = self[key] # force fun eval if relevant
            num = mfu.nw.identity(mfu.nw.shape(val)[0]) + 1.j*val
            denum = mfu.nw.identity(mfu.nw.shape(val)[0]) - 1.j*val
            newItem[key] = mfu.nw.dot(num, mfu.nw.invert(denum))
        return newItem
    def to_dTmat(self):
        raise NotImplementedError
    def to_dKmat(self):
        return self

    def to_dXSmat(self):
        raise NotImplementedError
    def to_dEPhaseMat(self):
        raise NotImplementedError
    def to_dUniOpMat(self):
        raise NotImplementedError

class dTmat(dMat):
    def __init__(self, d={}, asymCal=None, sourceStr=""):
        dMat.__init__(self, d, asymCal, sourceStr)
        self.chartTitle = "T matrix"

    def to_dSmat(self):
        raise NotImplementedError
    def to_dTmat(self):
        raise NotImplementedError
    def to_dKmat(self):
        raise NotImplementedError

    def to_dXSmat(self):
        raise NotImplementedError
    def to_dEPhaseMat(self):
        raise NotImplementedError
    def to_dUniOpMat(self):
        raise NotImplementedError

Smat = 0
Kmat = 1
Tmat = 2
def getDiscreteScatteringMatrix(matType, matDict, asymCal, sourceStr=""):
    if matType == Smat:
        return dSmat(matDict, asymCal, sourceStr)
    elif matType == Kmat:
        return dKmat(matDict, asymCal, sourceStr)
    elif matType == Tmat:
        return dTmat(matDict, asymCal, sourceStr)
    else:
        raise Exception("Non-recognised matrix type.")
