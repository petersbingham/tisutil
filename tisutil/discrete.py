import matfuncutil as mfu
from channelutil.units import *

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
        newItem = self._createNewItem(newUnits)
        self._initNewItem(newItem)
        for k,v in self.iteritems():
            newItem[k*fac] = v
        return newItem

class dVec(mfu.dVec, dBase):
    def _getReductionContainer(self):
        return dVal(units=self.units)

class dMat(mfu.dMat, dBase):
    def _getReductionContainer(self):
        return dVec(units=self.units)


class dSmat(dMat):
    def to_dTmat(self):
        newItem = self._createNewItem(self.units)
        self._initNewItem(newItem)
        for k,v in self.iteritems():
            newItem[k] = v - mfu.nw.identity(mfu.nw.shape(v)[0])
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
    def to_dTmat(self):
        raise NotImplementedError
    def to_dSmat(self):
        newItem = self._createNewItem(self.units)
        self._initNewItem(newItem)
        for k,v in self.iteritems():
            num = mfu.nw.identity(mfu.nw.shape(v)[0]) + 1.j*v
            denum = mfu.nw.identity(mfu.nw.shape(v)[0]) - 1.j*v
            newItem[k] = num / denum
        return newItem
    def to_dXSmat(self):
        raise NotImplementedError

class dTmat(dMat):
    def to_dSmat(self):
        raise NotImplementedError
    def to_dKmat(self):
        raise NotImplementedError
    def to_dXSmat(self):
        raise NotImplementedError


def usePythonTypes_d(dps=mfu.nw.dps_default_python):
    mfu.usePythonTypes(dps)

def useMpmathTypes_d(dps=mfu.nw.dps_default_mpmath):
    mfu.useMpmathTypes(dps)
