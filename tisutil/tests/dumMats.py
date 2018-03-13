from tisutil import discrete as dis

def rowOffsetColGain_zeroImag(sz=100):
    d = dis.dSmats()
    _rowOffsetColGain_zeroImag(d,sz)
    return d

def rowOffsetColGain_posImag(sz=100):
    d = dis.dmats()
    _rowOffsetColGain_posImag(d,sz)
    return d

def rowOffsetColGain_negImag(sz=100):
    d = dis.dmats()
    _rowOffsetColGain_negImag(d,sz)
    return d

def rowOffsetColGain_posNegImag(rg=10):
    d = dis.dSmats()
    _rowOffsetColGain_posNegImag(d,rg)
    return d

def rowOffsetColGain_zeroImag_Smat(sz=100):
    d = dis.dSmats()
    _rowOffsetColGain_zeroImag(d,sz)
    return d


def _rowOffsetColGain_posNegImag(d, rg):
    for i in range(-rg,rg+1):
        d[float(i)-float(i)*1j] = dis.mfu.nw.matrix([[float(i), 2*float(i)], 
                                                     [10+float(i), 10+2*float(i)]])
    return d

def _rowOffsetColGain_zeroImag(d, sz):
    for i in range(sz):
        d[float(i)] = dis.mfu.nw.matrix([[float(i), 2*float(i)], 
                                         [10+float(i), 10+2*float(i)]])
    return d

def _rowOffsetColGain_posImag(d, sz):
    for i in range(sz):
        d[float(i)+1j] = dis.mfu.nw.matrix([[float(i), 2*float(i)], 
                                            [10+float(i), 10+2*float(i)]])
    return d

def _rowOffsetColGain_negImag(d, sz):
    for i in range(sz):
        d[float(i)-1j] = dis.mfu.nw.matrix([[float(i), 2*float(i)], 
                                            [10+float(i), 10+2*float(i)]])
    return d
