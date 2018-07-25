import matfuncutil as mfu
from channelutil.units import *
import copy

class dBase:
    def convert_ene_units(self, new_x_units):
        if new_x_units == self.x_units:
            return self
        elif self.x_units==rydbergs:
            if new_x_units==hartrees:
                fac = rydbergs_to_hartrees
            elif new_x_units==eVs:
                fac = rydbergs_to_eVs
            else:
                raise Exception("Unknown conversion")
        elif self.x_units==hartrees:
            if new_x_units==rydbergs:
                fac = 1./rydbergs_to_hartrees
            elif new_x_units==eVs:
                fac = hartrees_to_eVs
            else:
                raise Exception("Unknown conversion")
        elif self.x_units==eVs:
            if new_x_units==rydbergs:
                fac = 1./rydbergs_to_eVs
            elif new_x_units==hartrees:
                fac = 1./hartrees_to_eVs
            else:
                raise Exception("Unknown conversion")
        else:
            raise Exception("Unknown conversion")
        new_item = self._create_new_item(x_units=new_x_units)
        self._init_new_item(new_item)
        for ene,val in self.iteritems():
            if not callable(val):
                new_item[ene*fac] = val
            else:
                new_item[ene*fac] = lambda new_ene: val(new_ene/fac)
        return new_item

class dSca(mfu.dSca, dBase):
    pass

class dXSsca(dSca):
    def __init__(self, d=None, x_units=None, source_str=None):
        y_units = "bohr$^2$"
        x_plotlbl = "Energy"
        y_plotlbl = "Cross Section"
        chart_title = "Cross Section"
        dSca.__init__(self, d, x_units, y_units, chart_title, x_plotlbl,
                      y_plotlbl, source_str)

class dEPhaseSca(dSca):
    def __init__(self, d=None, x_units=None, source_str=None):
        y_units = "radians"
        x_plotlbl = "Energy"
        y_plotlbl = "Eigenphase Sum"
        chart_title = "Eigenphase Sum"
        dSca.__init__(self, d, x_units, y_units, chart_title, x_plotlbl,
                      y_plotlbl, source_str)

class dVec(mfu.dVec, dBase):
    def _get_reduction_container(self):
        return dSca({}, self.x_units, self.y_units, self.chart_title,
                    self.x_plotlbl, self.y_plotlbl, self.source_str)

class dMat(mfu.dMat, dBase):
    def __init__(self, d=None, asymcalc=None, y_units=None, chart_title="",
                 x_plotlbl="", y_plotlbl="", source_str=""):
        mfu.dMat.__init__(self, d, 
                          None if asymcalc is None else asymcalc.get_units(), 
                          y_units, chart_title, x_plotlbl, y_plotlbl,
                          source_str)
        self.asymcalc = asymcalc
    def get_check_str(self):
        return mfu.dMat.get_check_str(self) + "\n" + str(self.asymcalc)
    def _get_reduction_container(self):
        return dVec({}, self.x_units, self.y_units, self.chart_title,
                    self.x_plotlbl, self.y_plotlbl, self.source_str)
    def _create_new_item(self, x_units=None, new_type=None):
        asymcalc = copy.deepcopy(self.asymcalc)
        if x_units is not None:
            asymcalc.units = x_units
        if new_type is None:
            new_type = type(self)
        new_item = new_type({}, asymcalc, self.source_str)
        return new_item
    def _I(self):
        return mfu.nw.identity(self.valshape()[0])
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
    def _get_dSca(self):
        return dSca(x_units=self.asymcalc.get_units(),
                    source_str=self.source_str)
    def _get_dVec(self):
        return dVec(x_units=self.asymcalc.get_units(),
                    source_str=self.source_str)
    def _get_dMat(self):
        return dMat(asymcalc=self.asymcalc, source_str=self.source_str)

class dSmat(dMat):
    def __init__(self, d=None, asymcalc=None, source_str=""):
        dMat.__init__(self, d, asymcalc, None, "S matrix", "Energy", "S matrix",
                      source_str)

    def to_dSmat(self):
        return self
    def to_dTmat(self):
        return self._convert(dTmat,
                             lambda val: val - self._I())
    def to_dKmat(self):
        return self._convert(dKmat,
                             lambda val: 1.j*self._I() - 1.j*val,
                             lambda val: self._I() + val)
    def to_dXSmat(self):
        return self.to_dTmat().to_dXSmat()
    def to_dXSsca(self):
        return self.to_dXSmat().to_dXSsca()

    def to_dEPhaseMat(self):
        return self.to_dKmat().to_dEPhaseMat()
    def to_dEPhaseSca(self):
        return self.to_dEPhaseMat().to_dEPhaseSca()

    def to_dUniOpMat(self):
        return self.unitary_op()

class dKmat(dMat):
    def __init__(self, d=None, asymcalc=None, source_str=""):
        dMat.__init__(self, d, asymcalc, None, "K matrix", "Energy", "K matrix",
                      source_str)

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
    def to_dXSsca(self):
        return self.to_dXSmat().to_dXSsca()

    def to_dEPhaseMat(self):
        new_item = self._create_new_item(new_type=dEPhasemat)
        self._init_new_item(new_item)
        dKmat_diag = self.diagonalise()
        for ene in dKmat_diag:
            val = dKmat_diag[ene] # force fun eval if relevant
            new_item[ene] = mfu.nw.arctan(val)
        return new_item
    def to_dEPhaseSca(self):
        return self.to_dEPhaseMat().to_dEPhaseSca()

class dTmat(dMat):
    def __init__(self, d=None, asymcalc=None, source_str=""):
        dMat.__init__(self, d, asymcalc, None, "T matrix", "Energy", "T matrix",
                      source_str)

    def to_dSmat(self):
        return self._convert(dSmat,
                             lambda val: val + self._I())
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
    def to_dXSsca(self):
        return self.to_dXSmat().to_dXSsca()

    def to_dEPhaseMat(self):
        return self.to_dKmat().to_dEPhaseMat()
    def to_dEPhaseSca(self):
        return self.to_dEPhaseMat().to_dEPhaseSca()

class dEPhasemat(dMat):
    def __init__(self, d=None, asymcalc=None, source_str=""):
        dMat.__init__(self, d, asymcalc, "radians", "Eigenphase Matrix",
                      "Energy", "Eigenphase", source_str)

    def to_dEPhaseSca(self):
        new_item = dEPhaseSca({}, self.x_units, self.source_str)
        self._init_new_item(new_item)
        for ene in self:
            val = self[ene] # force fun eval if relevant
            new_item[ene] = mfu.nw.sum_elements(val)
        return new_item

class dXSmat(dMat):
    def __init__(self, d=None, asymcalc=None, source_str=""):
        dMat.__init__(self, d, asymcalc, "bohr$^2$", "Cross Section", "Energy",
                      "Cross Section", source_str)

    def to_dXSsca(self):
        new_item = dXSsca({}, self.x_units, self.source_str)
        self._init_new_item(new_item)
        for ene in self:
            val = self[ene] # force fun eval if relevant
            new_item[ene] = mfu.nw.sum_elements(val)
        return new_item

class dFin(dMat):
    def __init__(self, d=None, asymcalc=None, source_str=""):
        dMat.__init__(self, d, asymcalc, None, "Fin", "Energy", "Fin",
                      source_str)

Smat = 0
Kmat = 1
Tmat = 2
def get_discrete_scattering_matrix(mat_type, mat_dict, asymcalc, source_str):
    if mat_type == Smat:
        return dSmat(mat_dict, asymcalc, source_str)
    elif mat_type == Kmat:
        return dKmat(mat_dict, asymcalc, source_str)
    elif mat_type == Tmat:
        return dTmat(mat_dict, asymcalc, source_str)
    else:
        raise Exception("Non-recognised matrix type.")
