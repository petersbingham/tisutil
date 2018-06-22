from discrete import *

class cMat(mfu.cMat):
    def __init__(self, fun_ref, asymcalc=None, source_str=""):
        mfu.cMat.__init__(self, fun_ref, 
                          None if asymcalc is None else asymcalc.get_units(),
                          source_str)
        self.asymcalc = asymcalc

class cSmat(cMat):
    def _get_discrete_container(self):
        return dSmat(asymcalc=self.asymcalc)

class cKmat(cMat):
    def _get_discrete_container(self):
        return dKmat(asymcalc=self.asymcalc)

class cTmat(cMat):
    def _get_discrete_container(self):
        return dTmat(asymcalc=self.asymcalc)

def get_continuous_scattering_matrix(mat_type, fun_ref, asymcalc, source_str=""):
    if mat_type == Smat:
        return cSmat(fun_ref, asymcalc, source_str)
    elif mat_type == Kmat:
        return cKmat(fun_ref, asymcalc, source_str)
    elif mat_type == Tmat:
        return cTmat(fun_ref, asymcalc, source_str)
    else:
        raise Exception("Non-recognised matrix type.")

# k as in wavenumber, not K-matrix
class cPolykmat(mfu.cSympyPolyMat):
    def __init__(self, sym_mat, sym_var, asymcalc, source_str=""):
        mfu.cMat.__init__(self,
          lambda ene: mfu.nw.from_sympy_matrix(sym_mat.subs(sym_var,
                                                            asymcalc.fk(ene))),
          asymcalc.get_units(), source_str)
        self.asymcalc = asymcalc
        self.sym_mat = sym_mat
        self.sym_var = sym_var

class cPolyFin(cPolykmat):
    def _get_discrete_container(self):
        return dFin(asymcalc=self.asymcalc)

class cPolySmat(cPolykmat):
    def _get_discrete_container(self):
        return dSmat(asymcalc=self.asymcalc)
