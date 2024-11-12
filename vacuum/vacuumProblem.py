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
    
    def solve(self, mpol: int, ntor: int, lrad: int):
        self.preset(mpol, ntor)

    def preset(self, mpol, ntor):
        self._set_resolution(mpol, ntor)
        
    def _set_resolution(self, mpol, ntor):
        import preset_fortran
        self.mpol, self.ntor = mpol, ntor
        self.mn = 1 + self.ntor + self.mpol*(2*self.ntor+1)
        self.xm, self.xn = preset_fortran.setresolution(
            self.mpol, self.ntor, self.nfp, self.mn
        )
        # extra-enhanced resolution for metrics
        self.mpol_metric, self.ntor_metric = 4*self.mpol, 4*self.ntor
        self.mn_metric = 1 + self.ntor_metric + self.mpol_metric*(2*self.ntor_metric+1)
        self.xm_metric, self.xn_metric = preset_fortran.setresolution(
            self.mpol_metric, self.ntor_metric, self.nfp, self.mn_metric
        )
        



if __name__ == '__main__':
    pass
