from discrete import *

def _set_chart_title_for_new(item, new_item):
    # Chart title has been supplemented in the client.
    # This is a bit of a hack since title set in discrete for tisutil.
    if item.chart_title != "":
        new_item.supplement_chart_title(item.chart_title)

class cMat(mfu.cMat):
    def __init__(self, fun_ref, asymcalc=None, source_str=""):
        mfu.cMat.__init__(self, fun_ref, 
                          None if asymcalc is None else asymcalc.get_units(),
                          source_str=source_str)
        self.asymcalc = asymcalc

class cSmat(cMat):
    def _get_discrete_container(self):
        dsmat = dSmat({}, self.asymcalc, self.source_str)
        _set_chart_title_for_new(self, dsmat)
        return dsmat

class cKmat(cMat):
    def _get_discrete_container(self):
        dkmat = dKmat({}, self.asymcalc, self.source_str)
        _set_chart_title_for_new(self, dkmat)
        return dkmat

class cTmat(cMat):
    def _get_discrete_container(self):
        dtmat = dTmat({}, self.asymcalc, self.source_str)
        _set_chart_title_for_new(self, dtmat)
        return dtmat

def get_continuous_scattering_matrix(mat_type, fun_ref, asymcalc, source_str):
    if mat_type == Smat:
        return cSmat(fun_ref, asymcalc, source_str)
    elif mat_type == Kmat:
        return cKmat(fun_ref, asymcalc, source_str)
    elif mat_type == Tmat:
        return cTmat(fun_ref, asymcalc, source_str)
    else:
        raise Exception("Non-recognised matrix type.")

# k as in wavenumber, not K-matrix
class cMatSympypolyk(mfu.cMatSympypoly):
    def __init__(self, sym_mat, sym_var, asymcalc, source_str=""):
        mfu.cMat.__init__(self,
          lambda ene: mfu.nw.from_sympy_matrix(sym_mat.subs(sym_var,
                                                            asymcalc.fk(ene))),
          asymcalc.get_units(), source_str=source_str)
        self.asymcalc = asymcalc
        self.sym_mat = sym_mat
        self.sym_var = sym_var

class cFinMatSympypolyk(cMatSympypolyk):
    def _get_discrete_container(self):
        dfin = dFin({}, self.asymcalc, self.source_str)
        _set_chart_title_for_new(self, dfin)
        return dfin

class cSMatSympypolyk(cMatSympypolyk):
    def _get_discrete_container(self):
        dsmat = dSmat({}, self.asymcalc, self.source_str)
        _set_chart_title_for_new(self, dsmat)
        return dsmat
