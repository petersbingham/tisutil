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

class dSca(mfu.dSca, dBase):
    pass

class dTotXSval(dSca):
    def __init__(self, d={}, units=None, source_str=""):
        dSca.__init__(self, d, units, source_str)
        self.chart_title = "Total Cross Section"

class dVec(mfu.dVec, dBase):
    def _get_reduction_container(self):
        return dSca(units=self.units)

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
    def _create_new_item(self, units=None, new_type=None):
        asymcalc = copy.deepcopy(self.asymcalc)
        if units is not None:
            asymcalc.units = units
        if new_type is None:
            new_type = type(self)
        new_item = new_type(asymcalc=asymcalc, source_str=self.source_str)
        return new_item
    def _I(self):
        return mfu.nw.identity(self._get_size())
    def _convert(self, new_type, num, denum=None):
        new_item = self._create_new_item(new_type=new_type)
        self._init_new_item(new_item)
        for ene in self:
            val = self[ene] # force fun eval if relevant
            if denum is not None:
                new_item[ene] = mfu.nw.dot(num(val), mfu.nw.invert(denum(val)))
            else:
                new_item[ene] = num(val)
        return new_item

class dSmat(dMat):
    def __init__(self, d={}, asymcalc=None, source_str=""):
        dMat.__init__(self, d, asymcalc, source_str)
        self.chart_title = "S matrix"

    def to_dSmat(self):
        return self
    def to_dTmat(self):
        return self._convert(dTmat,
                             lambda val: self._I() - val)
    def to_dKmat(self):
        return self._convert(dKmat,
                             lambda val: 1.j*self._I() - 1.j*val,
                             lambda val: self._I() + val)
    def to_dXSmat(self):
        return self.to_dTmat().to_dXSmat()
    def to_dEPhaseMat(self):
        raise NotImplementedError
    def to_dUniOpMat(self):
        return self.unitary_op()

class dKmat(dMat):
    def __init__(self, d={}, asymcalc=None, source_str=""):
        dMat.__init__(self, d, asymcalc, source_str)
        self.chart_title = "K matrix"

    def to_dSmat(self):
        return self._convert(dSmat,
                             lambda val: self._I() + 1.j*val,
                             lambda val: self._I() - 1.j*val)
    def to_dTmat(self):
        return self._convert(dTmat,
                             lambda val: 2.j*val,
                             lambda val: self._I() - 1.j*val)
    def to_dKmat(self):
        return self

    def to_dXSmat(self):
        return self.to_dTmat().to_dXSmat()
    def to_dEPhaseMat(self):
        raise NotImplementedError

class dTmat(dMat):
    def __init__(self, d={}, asymcalc=None, source_str=""):
        dMat.__init__(self, d, asymcalc, source_str)
        self.chart_title = "T matrix"

    def to_dSmat(self):
        return self._convert(dSmat,
                             lambda val: self._I() - val)
    def to_dTmat(self):
        return self
    def to_dKmat(self):
        return self._convert(dKmat,
                             lambda val: -1.j*val,
                             lambda val: 2*self._I() + val)

    def to_dXSmat(self):
        new_item = self._create_new_item(new_type=dXSmat)
        self._init_new_item(new_item)
        for ene in self:
            val = self[ene] # force fun eval if relevant
            def elementFunction(_, ch, el_val):
                k = self.asymcalc.k(ch, ene)
                a = mfu.nw.pi/mfu.nw.pow(k,2.)
                num = 2.*self.asymcalc.tot_spin() + 1.
                denum = 2.*(2.*self.asymcalc.targ_spins(ch) + 1.)
                b = num / denum
                c = mfu.nw.pow(mfu.nw.abs(el_val),2.)
                return a * b * c
            new_item[ene] = mfu.nw.apply_fun_to_elements(val, elementFunction)
        return new_item
    def to_dEPhaseMat(self):
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
