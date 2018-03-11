    def findRoot(self, start):
        def eneEqu(self, ene):
            self.setEnergy(ene, False)
            return self._denum(False)
        return mpm.findroot(lambda ene: eneEqu(self, ene), start)
    
    def plotPoles(self, Rs, Is, real):
        self.i = 0
        @np.vectorize
        def eneEqu(self, R, I):
            print self.i
            self.i += 1
            self.setEnergy(R+I*1.0j, False)
            return 1.0 / self._denum(False)
        plot3D.plot(Rs, Is, lambda R, I: eneEqu(self, R, I), real, "Energy", "1/Denum")
        
    # Matrix that has absolute values of elements.
    def abs(self):
        return Smat.AbsSmat(self)
    
    class UniOpSmat(sm.mat):
        def __init__(self, Smat):
            sm.mat.__init__(self, 2, num.PRECISION)
            self.Smat = Smat
            self.calculate()
        def calculate(self):
            m1 = self.Smat.getMatrix()
            m2 = m1.transpose().conjugate()
            self.uniOpS = m1 * m2
        def _getRow(self, i):
            return [self.uniOpS[i,0], self.uniOpS[i,1]]
      
    def isUnitary(self):
        if not gu.complexCompare(self.uniOpSmat[0][0],1.0):
            return False
        if not gu.complexCompare(self.uniOpSmat[0][1],0.0):
            return False
        if not gu.complexCompare(self.uniOpSmat[1][0],0.0):
            return False
        if not gu.complexCompare(self.uniOpSmat[1][1],1.0):
            return False
        return True