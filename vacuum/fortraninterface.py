#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# fortraninterface.py


import vacuum_fortran
import numpy as np
from typing import Tuple


# basefunction.f90
def get_zernike(r: float, lrad: int, mpol: int) -> np.ndarray:
    return vacuum_fortran.get_zernike(r, lrad, mpol)


def get_zernike_d2(r: float, lrad: int, mpol: int) -> np.ndarray:
    return vacuum_fortran.get_zernike_d2(r, lrad, mpol)


def get_zernike_rm(r: float, lrad: int, mpol: int) -> np.ndarray:
    return vacuum_fortran.get_zernike_rm(r, lrad, mpol)



# preset.f90
def get_resolution(mpol: int, ntor: int, nfp: int) -> Tuple[np.ndarray]:
    mn = 1 + ntor + mpol*(2*ntor+1)
    return vacuum_fortran.get_resolution(mpol, ntor, nfp, mn)

def get_NAdof(mpol: int, ntor: int, lrad: int, stellsym: bool) -> int:
    mn = 1 + ntor + mpol*(2*ntor+1)
    return vacuum_fortran.get_NAdof(mpol, ntor, mn, lrad, int(stellsym))