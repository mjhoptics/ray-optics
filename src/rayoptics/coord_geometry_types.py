#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" type hints for for vectors and matrices

These type hints are provided to ensure a consist convention of distinguishing between numpy arrays and array-like arrays.

2D, 3D and 4D cases are provided.

VecNd is used for coordinates
DirNd is used for vector directions, unit length
MatNd is a N x N matrix
TfmNd is used to package together a rotation matrix and translation vector

Ray3d, a tuple of a Vec3d and Dir3d is defined by a point on the ray and 
the ray direction.

Created on Mon Jun  2 10:48:54 2025

.. codeauthor: Michael J. Hayford
"""
import numpy as np
import numpy.typing as npt

Vec2d = npt.NDArray[np.double]
Dir2d = npt.NDArray[np.double]
Mat2d = npt.NDArray[np.double]
Tfm2d = tuple[Mat2d|None, Vec2d]

Vec3d = npt.NDArray[np.double]
Dir3d = npt.NDArray[np.double]
Mat3d = npt.NDArray[np.double]
Tfm3d = tuple[Mat3d|None, Vec3d]
Ray3d = tuple[Vec3d, Dir3d]

Vec4d = npt.NDArray[np.double]
Dir4d = npt.NDArray[np.double]
Mat4d = npt.NDArray[np.double]
Tfm4d = tuple[Mat4d|None, Vec4d]

V2d = npt.ArrayLike
M2d = npt.ArrayLike
T2d = tuple[M2d|None, V2d]

V3d = npt.ArrayLike
D3d = npt.ArrayLike
M3d = npt.ArrayLike
T3d = tuple[M3d|None, V3d]

V4d = npt.ArrayLike
M4d = npt.ArrayLike
T4d = tuple[M4d|None, V4d]
