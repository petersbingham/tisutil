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
        newItem = self._createNewItem(units=newUnits)
        self._initNewItem(newItem)
        for ene,val in self.iteritems():
            newItem[ene*fac] = val
        return newItem

class dVal(mfu.dVal, dBase):
    pass

class dTotXSval(dVal):
    def __init__(self, d={}, units=None, sourceStr=""):
        dVal.__init__(self, d, units, sourceStr)
        self.chartTitle = "Total Cross Section"

class dVec(mfu.dVec, dBase):
    def _getReductionContainer(self):
        return dVal(units=self.units)

class dMat(mfu.dMat, dBase):
    def __init__(self, d={}, asymCal=None, sourceStr=""):
        mfu.dMat.__init__(self, d, 
                          None if asymCal is None else asymCal.getUnits(), 
                          sourceStr)
        self.asymCal = asymCal
    def getCheckStr(self):
        return mfu.dMat.getCheckStr(self) + "\n" + str(self.asymCal)
    def _getReductionContainer(self):
        return dVec(units=self.units)
    def _createNewItem(self, units=None, newType=None):
        asymCal = copy.deepcopy(self.asymCal)
        if units is not None:
            asymCal.units = units
        if newType is None:
            newType = type(self)
        newItem = newType(asymCal=asymCal, sourceStr=self.sourceStr)
        return newItem

class dSmat(dMat):
    def __init__(self, d={}, asymCal=None, sourceStr=""):
        dMat.__init__(self, d, asymCal, sourceStr)
        self.chartTitle = "S matrix"

    def to_dSmat(self):
        return self
    def to_dTmat(self):
        newItem = self._createNewItem(newType=dTmat)
        self._initNewItem(newItem)
        for ene in self:
            val = self[ene] # force fun eval if relevant
            newItem[ene] = mfu.nw.identity(mfu.nw.shape(val)[0]) - val
        return newItem
    def to_dKmat(self):
        raise NotImplementedError

    def to_dXSmat(self):
        return self.to_dTmat().to_dXSmat()
    def to_dEPhaseMat(self):
        raise NotImplementedError
    def to_dUniOpMat(self):
        raise NotImplementedError

class dKmat(dMat):
    def __init__(self, d={}, asymCal=None, sourceStr=""):
        dMat.__init__(self, d, asymCal, sourceStr)
        self.chartTitle = "K matrix"

    def to_dSmat(self):
        newItem = self._createNewItem(newType=dSmat)
        self._initNewItem(newItem)
        for ene in self:
            val = self[ene] # force fun eval if relevant
            num = mfu.nw.identity(mfu.nw.shape(val)[0]) + 1.j*val
            denum = mfu.nw.identity(mfu.nw.shape(val)[0]) - 1.j*val
            newItem[ene] = mfu.nw.dot(num, mfu.nw.invert(denum))
        return newItem
    def to_dTmat(self):
        newItem = self._createNewItem(newType=dTmat)
        self._initNewItem(newItem)
        for ene in self:
            val = self[ene] # force fun eval if relevant
            num = 2.j*val
            denum = mfu.nw.identity(mfu.nw.shape(val)[0]) - 1.j*val
            newItem[ene] = mfu.nw.dot(num, mfu.nw.invert(denum))
        return newItem
    def to_dKmat(self):
        return self

    def to_dXSmat(self):
        return self.to_dTmat().to_dXSmat()
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
        return self
    def to_dKmat(self):
        raise NotImplementedError

    def to_dXSmat(self):
        newItem = self._createNewItem(newType=dXSmat)
        self._initNewItem(newItem)
        for ene in self:
            val = self[ene] # force fun eval if relevant
            def elementFunction(_, j, elVal):
                k = self.asymCal.k(j, ene)
                a = mfu.nw.pi/mfu.nw.pow(k,2.)
                b = mfu.nw.pow(mfu.nw.abs(elVal),2.)
                return a * b
            newItem[ene] = mfu.nw.applyFunToElements(val, elementFunction)
        return newItem
    def to_dEPhaseMat(self):
        raise NotImplementedError
    def to_dUniOpMat(self):
        raise NotImplementedError

class dXSmat(dMat):
    def __init__(self, d={}, asymCal=None, sourceStr=""):
        dMat.__init__(self, d, asymCal, sourceStr)
        self.chartTitle = "Cross Section"
    def to_dTotXSval(self):
        newItem = dTotXSval(units=self.units, sourceStr=self.sourceStr)
        self._initNewItem(newItem)
        for ene in self:
            val = self[ene] # force fun eval if relevant
            newItem[ene] = mfu.nw.sumElements(val)
        return newItem

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
