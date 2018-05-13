import matfuncutil as mfu
from channelutil.units import *
import copy

class dBase:
    def convert_units(self, new_units):
        if new_units == self.units:
            return self
        elif self.units==rydbergs:
            if new_units==hartrees:
                fac = rydbergs_to_hartrees
            elif new_units==eVs:
                fac = rydbergs_to_eVs
            else:
                raise Exception("Unknown conversion")
        elif self.units==hartrees:
            if new_units==rydbergs:
                fac = 1./rydbergs_to_hartrees
            elif new_units==eVs:
                fac = hartrees_to_eVs
            else:
                raise Exception("Unknown conversion")
        elif self.units==eVs:
            if new_units==rydbergs:
                fac = 1./rydbergs_to_eVs
            elif new_units==hartrees:
                fac = 1./hartrees_to_eVs
            else:
                raise Exception("Unknown conversion")
        else:
            raise Exception("Unknown conversion")
        new_item = self._create_new_item(units=new_units)
        self._init_new_item(new_item)
        for ene,val in self.iteritems():
            new_item[ene*fac] = val
        return new_item

class dVal(mfu.dVal, dBase):
    pass

class dTotXSval(dVal):
    def __init__(self, d={}, units=None, source_str=""):
        dVal.__init__(self, d, units, source_str)
        self.chart_title = "Total Cross Section"

class dVec(mfu.dVec, dBase):
    def _get_reduction_container(self):
        return dVal(units=self.units)

class dMat(mfu.dMat, dBase):
    def __init__(self, d={}, asymcalc=None, source_str=""):
        mfu.dMat.__init__(self, d, 
                          None if asymcalc is None else asymcalc.get_units(), 
                          source_str)
        self.asymcalc = asymcalc
    def get_check_str(self):
        return mfu.dMat.get_check_str(self) + "\n" + str(self.asymcalc)
    def _get_reduction_container(self):
        return dVec(units=self.units)
    def _create_new_item(self, units=None, newType=None):
        asymcalc = copy.deepcopy(self.asymcalc)
        if units is not None:
            asymcalc.units = units
        if newType is None:
            newType = type(self)
        new_item = newType(asymcalc=asymcalc, source_str=self.source_str)
        return new_item

class dSmat(dMat):
    def __init__(self, d={}, asymcalc=None, source_str=""):
        dMat.__init__(self, d, asymcalc, source_str)
        self.chart_title = "S matrix"

    def to_dSmat(self):
        return self
    def to_dTmat(self):
        new_item = self._create_new_item(newType=dTmat)
        self._init_new_item(new_item)
        for ene in self:
            val = self[ene] # force fun eval if relevant
            new_item[ene] = mfu.nw.identity(mfu.nw.shape(val)[0]) - val
        return new_item
    def to_dKmat(self):
        raise NotImplementedError

    def to_dXSmat(self):
        return self.to_dTmat().to_dXSmat()
    def to_dEPhaseMat(self):
        raise NotImplementedError
    def to_dUniOpMat(self):
        raise NotImplementedError

class dKmat(dMat):
    def __init__(self, d={}, asymcalc=None, source_str=""):
        dMat.__init__(self, d, asymcalc, source_str)
        self.chart_title = "K matrix"

    def to_dSmat(self):
        new_item = self._create_new_item(newType=dSmat)
        self._init_new_item(new_item)
        for ene in self:
            val = self[ene] # force fun eval if relevant
            num = mfu.nw.identity(mfu.nw.shape(val)[0]) + 1.j*val
            denum = mfu.nw.identity(mfu.nw.shape(val)[0]) - 1.j*val
            new_item[ene] = mfu.nw.dot(num, mfu.nw.invert(denum))
        return new_item
    def to_dTmat(self):
        new_item = self._create_new_item(newType=dTmat)
        self._init_new_item(new_item)
        for ene in self:
            val = self[ene] # force fun eval if relevant
            num = 2.j*val
            denum = mfu.nw.identity(mfu.nw.shape(val)[0]) - 1.j*val
            new_item[ene] = mfu.nw.dot(num, mfu.nw.invert(denum))
        return new_item
    def to_dKmat(self):
        return self

    def to_dXSmat(self):
        return self.to_dTmat().to_dXSmat()
    def to_dEPhaseMat(self):
        raise NotImplementedError
    def to_dUniOpMat(self):
        raise NotImplementedError

class dTmat(dMat):
    def __init__(self, d={}, asymcalc=None, source_str=""):
        dMat.__init__(self, d, asymcalc, source_str)
        self.chart_title = "T matrix"

    def to_dSmat(self):
        raise NotImplementedError
    def to_dTmat(self):
        return self
    def to_dKmat(self):
        raise NotImplementedError

    def to_dXSmat(self):
        new_item = self._create_new_item(newType=dXSmat)
        self._init_new_item(new_item)
        for ene in self:
            val = self[ene] # force fun eval if relevant
            def elementFunction(_, j, elVal):
                k = self.asymcalc.k(j, ene)
                a = mfu.nw.pi/mfu.nw.pow(k,2.)
                b = mfu.nw.pow(mfu.nw.abs(elVal),2.)
                return a * b
            new_item[ene] = mfu.nw.apply_fun_to_elements(val, elementFunction)
        return new_item
    def to_dEPhaseMat(self):
        raise NotImplementedError
    def to_dUniOpMat(self):
        raise NotImplementedError

class dXSmat(dMat):
    def __init__(self, d={}, asymcalc=None, source_str=""):
        dMat.__init__(self, d, asymcalc, source_str)
        self.chart_title = "Cross Section"
    def to_dTotXSval(self):
        new_item = dTotXSval(units=self.units, source_str=self.source_str)
        self._init_new_item(new_item)
        for ene in self:
            val = self[ene] # force fun eval if relevant
            new_item[ene] = mfu.nw.sum_elements(val)
        return new_item

Smat = 0
Kmat = 1
Tmat = 2
def get_discrete_scattering_matrix(mat_type, mat_dict, asymcalc, source_str=""):
    if mat_type == Smat:
        return dSmat(mat_dict, asymcalc, source_str)
    elif mat_type == Kmat:
        return dKmat(mat_dict, asymcalc, source_str)
    elif mat_type == Tmat:
        return dTmat(mat_dict, asymcalc, source_str)
    else:
        raise Exception("Non-recognised matrix type.")
