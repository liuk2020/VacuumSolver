#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# axis.py


import numpy as np
from typing import Tuple


class Axis:
    
    def __init__(self, nfp: int, ntor: int, rRe: np.ndarray, zIm: np.ndarray, rIm: np.ndarray, zRe: np.ndarray, stellSym: bool) -> None:
        self._nfp = nfp
        self.ntor = ntor
        self._stellSym = stellSym
        self._rRe = rRe
        self._zIm = zIm
        self._rIm = rIm
        self._zRe = zRe
        
    @property
    def nfp(self):
        return self._nfp
    
    @property
    def xn(self):
        return np.arange(self.ntor+1)
                              
    @property
    def stellSym(self):
        return self._stellSym
    
    @property
    def rRe(self):
        return self._rRe
    
    @property
    def zIm(self):
        return self._zIm
    
    @property
    def rIm(self):
        if self.stellSym:
            return np.zeros(self.ntor+1)
        else:
            return self._rIm
    
    @property
    def zRe(self):
        if self.stellSym:
            return np.zeros(self.ntor+1)
        else:
            return self._zRe
    
    def getRZ(self, zetaArr: np.ndarray) -> Tuple[np.ndarray]:
        if not isinstance(zetaArr, np.ndarray):
            try:
                zetaArr = np.array(zetaArr)
            except:
                zetaArr = np.array([zetaArr])
        angleMat = -self.nfp * np.dot(self.xn.reshape(-1,1), zetaArr.reshape(1,-1))
        rArr = 2 * (
            np.dot(self.rRe.reshape(1,-1), np.cos(angleMat)) - 
            np.dot(self.rIm.reshape(1,-1), np.sin(angleMat))
        )
        zArr = 2 * (
            np.dot(self.zRe.reshape(1,-1), np.cos(angleMat)) - 
            np.dot(self.zIm.reshape(1,-1), np.sin(angleMat))
        )
        rArr -= self.rRe[0]
        zArr -= self.zRe[0]
        try:
            m, n = zetaArr.shape
            return rArr.reshape(m, n), zArr.reshape(m, n)
        except:
            if isinstance(zetaArr, np.ndarray) and zetaArr.shape[0] == 1: 
                return rArr.flatten(), zArr.flatten()
            elif isinstance(zetaArr, np.ndarray) and zetaArr.shape[1] == 1: 
                return rArr.flatten(), zArr.flatten()
            else:
                return rArr, zArr


if __name__ == '__main__':
    pass