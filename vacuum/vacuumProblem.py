#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# vacuumProblem.py


import numpy as np
from .axis import Axis
from tfpy.geometry import Surface_cylindricalAngle


class Vacuum():
    
    def __init__(self, boundary: Surface_cylindricalAngle, axis: Axis=None) -> None:
        """
        Initialize the toroidal vacuum field!

        Parameters
        ----------
        boundary: Surface_cylindricalAngle
            The boundary of the torodial vacuum field
        """
        self.boundary = boundary
        if isinstance(axis, Axis):
            assert self.stellSym == axis.stellSym
            assert self.nfp == axis.nfp
            self.axis = axis
        else:
            self._init_axis()
        
    def _init_axis(self):
        self.axis(
            nfp = self.nfp,
            ntor = self.boundary.r.ntor,
            rRe = self.boundary.r._reArr[np.where(self.boundary.r.xm==0)[0]],
            zIm = self.boundary.z._imArr[np.where(self.boundary.z.xm==0)[0]],
            rIm = self.boundary.r._imArr[np.where(self.boundary.r.xm==0)[0]],
            zRe = self.boundary.z._reArr[np.where(self.boundary.z.xm==0)[0]],
            stellSym = self.stellSym
        )
    
    @property
    def stellSym(self) -> bool:
        return self.boundary
    
    @property
    def nfp(self):
        return self.boundary.nfp


if __name__ == '__main__':
    pass
