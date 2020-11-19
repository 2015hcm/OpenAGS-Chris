# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 11:43:52 2020

@author: chris
"""

from numpy import NaN, Inf, arange, asarray, array
import numpy as np
import xylib
from abc import ABC, abstractmethod

class Peak(ABC):
    def __init__(self):
        super().__init__()
        
    @abstractmethod
    def get_area(self):
        pass
    
    @abstractmethod
    def get_num_params(self):
        pass
    
    @abstractmethod
    def update_params(self, newParams):
        pass
    
    @abstractmethod
    def get_ydata(self, xdata):
        pass

class Background(ABC):
    def __init__(self):
        super().__init__()
        
    @abstractmethod
    def get_num_params(self):
        pass
    
    @abstractmethod
    def update_params(self, newParams):
        pass
    
    @abstractmethod
    def get_ydata(self, xdata):
        pass

class LinearBackground(Background):
    def __init__(self, slope=None, intercept = None, pointA = None, pointB = None):
        super().__init__()
        if slope != None and intercept != None:
            self.slope = slope
            self.intergept = intercept
        elif pointA != None and pointB != None:
            x1 = pointA[0]
            y1 = pointA[1]
            x2 = pointB[0]
            y2 = pointB[1]
            self.slope = (y2-y1) / (x2 - x1)
            self.intercept = y1 - self.slope * x1
    def get_num_params(self):
        return 2
    def update_params(self, newParams):
        self.slope = newParams[0]
        self.intercept = newParams[1]
    def get_ydata(self, xdata):
        return self.slope * xdata + self.intercept
    
def peakdet(v, delta, x = None):
    """
    Peak finder
    Converted from MATLAB script at http://billauer.co.il/peakdet.html
    
    Returns two arrays
    
    function [maxtab, mintab]=peakdet(v, delta, x)
    %PEAKDET Detect peaks in a vector
    %        [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
    %        maxima and minima ("peaks") in the vector V.
    %        MAXTAB and MINTAB consists of two columns. Column 1
    %        contains indices in V, and column 2 the found values.
    %      
    %        With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
    %        in MAXTAB and MINTAB are replaced with the corresponding
    %        X-values.
    %
    %        A point is considered a maximum peak if it has the maximal
    %        value, and was preceded (to the left) by a value lower by
    %        DELTA.
    
    % Eli Billauer, 3.4.05 (Explicitly not copyrighted).
    % This function is released to the public domain; Any use is allowed.
    
    """
    maxtab = []
    mintab = []
       
    if x is None:
        x = arange(len(v))
    
    v = asarray(v)
    
    mn, mx = Inf, -Inf
    mnpos, mxpos = NaN, NaN
    
    lookformax = True
    
    for i in arange(len(v)):
        this = v[i]
        if this > mx:
            mx = this
            mxpos = x[i]
        if this < mn:
            mn = this
            mnpos = x[i]
        
        if lookformax:
            if this < mx-delta:
                maxtab.append((mxpos, mx))
                mn = this
                mnpos = x[i]
                lookformax = False
        else:
            if this > mn+delta:
                mintab.append((mnpos, mn))
                mx = this
                mxpos = x[i]
                lookformax = True

    return maxtab

def var_mul(e1, var1, e2, var2): #assume these are independent variables, and multiply their variances together
    return var1*var2 + e1**2*var2 + e2**2 *var1

class SpectrumParser:
    """Parser for SPE files"""
    def __init__(self, fname):
        self.fname = fname
        if self.fname.split(".")[-1].lower() == "spe":
            self.spectrumfile = open(fname)
            self.speFile = True
        else:
            self.speFile = False
    def getValues(self):
        if self.speFile:
            text = self.spectrumfile.read()
            l = text.split('$')[1:]
            sections = dict()
            for section in l:
                sp = section.split('\n',1)
                sections[sp[0]] = sp[1]
            livetime, realtime = sections["MEAS_TIM:"].strip().split(" ")
            livetime = float(livetime)
            realtime = float(realtime)
            startind, endind = sections["DATA:"].split("\n")[0].split(" ")
            data = sections["DATA:"].split("\n")[1:-1]
            data = [int(i) for i in data]
            intercept, slope = sections["ENER_FIT:"].strip().split(" ")[:2]
            intercept = float(intercept)
            slope = float(slope)
            energies = [intercept + i*slope for i in range(int(startind),int(endind)+1)]
            cps = [total/livetime for total in data]
            return (livetime, realtime, np.array(energies), np.array(cps))
        else:
            d = xylib.load_file(self.fname, '')
            block = d.get_block(0)
            meta = block.meta
            metaDict = {}
            for i in range(meta.size()):
                key = meta.get_key(i)
                metaDict[key] = meta.get(key)
            ncol = block.get_column_count()
            nrow = block.get_point_count()
            data = [[block.get_column(i).get_value(j) for j in range(nrow)] for i in range(1, ncol+1)]
            livetime = metaDict["live time (s)"]
            realtime = metaDict["real time (s)"]
            return (livetime, realtime, np.array(data[0]), np.array(data[1])/float(livetime))
    def close(self):
        self.spectrumfile.close()
"""Functions to search for values in my data, used as utilities in many places. Modified binary search algorithm."""
def binary_search_find_nearest(l, e):
    upperBound = len(l)
    lowerBound = 0
    guess = (upperBound + lowerBound)//2
    while not (e < l[guess+1] and e > l[guess-1]):
        if e > l[guess]:
            lowerBound = guess + 1
        else:
            upperBound = guess - 1
        guess = (upperBound + lowerBound)//2
        if guess <= 2 or guess >= len(l)-2:
            break
    if e > l[guess]:
        guess += 1
    return guess

def binary_search_buried(l, e, i):
    upperBound = len(l)
    lowerBound = 0
    guess = (upperBound + lowerBound)//2
    while not (e < l[guess+1][i] and e > l[guess-1][i]):
        if e > l[guess][i]:
            lowerBound = guess + 1
        else:
            upperBound = guess - 1
        guess = (upperBound + lowerBound)//2
        if guess <= 2 or guess >= len(l)-2:
            break
    if e > l[guess][i]:
        guess += 1
    return guess
