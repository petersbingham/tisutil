import matfuncutil as mfu
import conversions as con

RYDs = 0
eVs = 1

class dBase:
    def convertUnits(self, newUnits):
        if newUnits == self.units:
            return self
        elif self.units==RYDs or self.units==eVs:
            if newUnits==RYDs:
                fac = 1./con.RYD_to_EV
            elif newUnits==eVs:
                fac = con.RYD_to_EV
            else:
                raise Exception("Unknown conversion")
        else:
            raise Exception("Unknown conversion")
        newItem = self._createNewItem(newUnits)
        self._initNewItem(newItem)
        for k,v in self.iteritems():
            newItem[k*fac] = v
        return newItem

class dvals(mfu.dvals, dBase):
    def __init__(self, d={}, units=RYDs):
        mfu.dvals.__init__(self, d, units)

class dvecs(mfu.dvecs, dBase):
    def __init__(self, d={}, units=RYDs):
        mfu.dvecs.__init__(self, d, units)

    def _getReductionContainer(self):
        return dvals(units=self.units)

class dmats(mfu.dmats, dBase):
    def __init__(self, d={}, units=RYDs):
        mfu.dmats.__init__(self, d, units)

    def _getReductionContainer(self):
        return dvecs(units=self.units)


class dSmats(dmats):
    def to_dTMats(self):
        newItem = self._createNewItem(self.units)
        self._initNewItem(newItem)
        for k,v in self.iteritems():
            newItem[k] = v - mfu.nw.identity(mfu.nw.shape(v)[0])
        return newItem
    def to_dKMats(self):
        raise NotImplementedError
    def to_dXSMats(self):
        raise NotImplementedError
    def to_dEPhaseMats(self):
        raise NotImplementedError
    def to_dUniOpMats(self):
        raise NotImplementedError

class dKmats(dmats):
    def to_dTMats(self):
        raise NotImplementedError
    def to_dSMats(self):
        raise NotImplementedError
    def to_dXSMats(self):
        raise NotImplementedError

class dTmats(dmats):
    def to_dSMats(self):
        raise NotImplementedError
    def to_dKMats(self):
        raise NotImplementedError
    def to_dXSMats(self):
        raise NotImplementedError
