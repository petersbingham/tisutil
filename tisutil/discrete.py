import numpy as np
import matplotlib.pyplot as plt
from cycler import cycler
import random

import pynumwrap as nw

import conversions as con

RYDs = 0
eVs = 1

class Base(dict):
    def __init__(self, d={}, units=RYDs):
        dict.__init__(self, d)
        self.units = units
        
        self.chartTitle = ""
        self.colourCycle = ['red', 'green', 'blue', 'purple']
        self.legPrefix = ""
        self.useMarker = False
        self.xsize = None
        self.ysize = None
        
        self.sigFigs = 6

    
    def convertUnits(self, units):
        if units == self.units:
            return self
        elif self.units==RYDs or self.units==eVs:
            if units==RYDs:
                fac = 1./con.RYD_to_EV
            elif units==eVs:
                fac = con.RYD_to_EV
            else:
                raise("Unknown conversion")
        else:
            raise("Unknown conversion")
        newItem = self._createNewItem(units)
        self._initNewItem(newItem)
        for k,v in self.iteritems():
            newItem[k*fac] = v
        return newItem

class Smats(mats):
    def toTMats(self):
        newItem = self._createNewItem(self.units)
        self._initNewItem(newItem)
        for k,v in self.iteritems():
            newItem[k] = v - nw.identity(nw.shape(v)[0])
        return newItem
    def toKMats(self):
        raise NotImplementedError
    def toXSMats(self):
        raise NotImplementedError
    def toEPhaseMats(self):
        raise NotImplementedError
    def toUniOpMats(self):
        raise NotImplementedError

class Kmats(mats):
    def toTMats(self):
        raise NotImplementedError
    def toSMats(self):
        raise NotImplementedError
    def toXSMats(self):
        raise NotImplementedError

class Tmats(mats):
    def toSMats(self):
        raise NotImplementedError
    def toKMats(self):
        raise NotImplementedError
    def toXSMats(self):
        raise NotImplementedError
