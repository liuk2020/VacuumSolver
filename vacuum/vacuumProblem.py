#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# vacuumProblem.py


from tfpy.geometry import Surface_cylindricalAngle


class Vacuum():
    
    def __init__(self, boundary: Surface_cylindricalAngle) -> None:
        """
        Initialize the toroidal vacuum field!

        Parameters
        ----------
        boundary: Surface_cylindricalAngle
            The boundary of the torodial vacuum field
        """
        self.boundary = boundary


if __name__ == '__main__':
    pass
