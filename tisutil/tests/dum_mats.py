import tisutil as tu
import channelutil as cu

def row_offset_col_gain_zeroimag(sz=100):
    d = tu.dSmat()
    _row_offset_col_gain_zeroimag(d,sz)
    return d

def row_offset_col_gain_posImag(sz=100):
    d = tu.dMat()
    _row_offset_col_gain_posImag(d,sz)
    return d

def row_offset_col_gain_negImag(sz=100):
    d = tu.dMat()
    _row_offset_col_gain_negImag(d,sz)
    return d

def row_offset_col_gain_posNegImag(rg=10):
    d = tu.dSmat(asymcalc=cu.AsymCalc(tu.rydbergs,[0,0]))
    _row_offset_col_gain_posNegImag(d,rg)
    return d

def row_offset_col_gain_zeroimag_Smat(sz=100):
    d = tu.dSmat()
    _row_offset_col_gain_zeroimag(d,sz)
    return d


def _row_offset_col_gain_posNegImag(d, rg):
    for i in range(-rg,rg+1):
        d[float(i)-float(i)*1j] = tu.mfu.nw.matrix([[float(i), 2*float(i)], 
                                                     [10+float(i), 10+2*float(i)]])
    return d

def _row_offset_col_gain_zeroimag(d, sz):
    for i in range(sz):
        d[float(i)] = tu.mfu.nw.matrix([[float(i), 2*float(i)], 
                                         [10+float(i), 10+2*float(i)]])
    return d

def _row_offset_col_gain_posImag(d, sz):
    for i in range(sz):
        d[float(i)+1j] = tu.mfu.nw.matrix([[float(i), 2*float(i)], 
                                            [10+float(i), 10+2*float(i)]])
    return d

def _row_offset_col_gain_negImag(d, sz):
    for i in range(sz):
        d[float(i)-1j] = tu.mfu.nw.matrix([[float(i), 2*float(i)], 
                                            [10+float(i), 10+2*float(i)]])
    return d
