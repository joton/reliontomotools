#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 28 10:23:52 2021

@author: joton
"""

import numpy as np
# from transforms3d.euler import euler2mat, mat2euler


# % Euler transformations


def rotX(angle):

    c = np.cos(np.deg2rad(angle))
    s = np.sin(np.deg2rad(angle))

    A = np.array([[1, 0, 0, 0],
                  [0, c,-s, 0],
                  [0, s, c, 0],
                  [0, 0, 0, 1]])
    return A


def rotY(angle):

    c = np.cos(np.deg2rad(angle))
    s = np.sin(np.deg2rad(angle))

    A = np.array([[c, 0, s, 0],
                  [0, 1, 0, 0],
                  [-s,0, c, 0],
                  [0, 0, 0, 1]])
    return A


def rotZ(angle):

    c = np.cos(np.deg2rad(angle))
    s = np.sin(np.deg2rad(angle))

    A = np.array([[c,-s, 0, 0],
                  [s, c, 0, 0],
                  [0, 0, 1, 0],
                  [0, 0, 0, 1.]])
    return A


def getShiftMatrix(shifts):

    if not isinstance(shifts, np.ndarray):
        shifts = np.array(shifts)

    shiftM = np.identity(4)
    shiftM[0:3, 3] = shifts
    return shiftM


# def getTransMatrix(angles=np.array([0, 0, 0]), shifts=np.array([0, 0, 0])):
#     """
#     Parameters
#     ----------
#     angles : 1D array
#         with Rot, Tilt and Psi Relion angles
#     shifts : 1D array
#         Xshift, Yshift and Zshift

#     Returns
#     -------
#     Trans matrix  : 2D array with transformation matrix which applies from
#                     reference to particle
#     """

#     if not isinstance(angles, np.ndarray):
#         angles = np.array(angles)
#     if not isinstance(shifts, np.ndarray):
#         shifts = np.array(shifts)

#     transMat = np.zeros((4, 4))
#     transMat[3, 3] = 1
#     transMat[:3, :3] = euler2mat(*(np.pi/180*angles), 'sxyz')
#     transMat[:3, 3] = shifts

#     return transMat


# def getTransFromMatrix(transMatrix):

#     shifts = transMatrix[0:3, 3]

#     angles = np.array(mat2euler(transMatrix[0:3, 0:3], 'sxyz'))*180/np.pi

#     if angles[1] < 0:
#         angles[[0, 2]] += 180
#         angles[1] *= -1

#     return angles, shifts

