# -*- coding: utf-8 -*-
"""
Created on Sun Oct  4 19:01:50 2020

@author: malte
"""

import numpy as np
# from scipy.stats import special_ortho_group
import sys
import math as m


def symToMatrix(symTensor):
    # xx,xy,xz,yy,yz,zz
    Matrix = np.array([[symTensor[0], symTensor[1], symTensor[2]],
                       [symTensor[1], symTensor[3], symTensor[4]],
                       [symTensor[2], symTensor[4], symTensor[5]]])
    return Matrix


def symToMatrixPVPoly(symTensor):
    # xx,xy,xz,yy,yz,zz
    Matrix = np.array([[symTensor[0], symTensor[3], symTensor[4]],
                       [symTensor[3], symTensor[1], symTensor[5]],
                       [symTensor[4], symTensor[5], symTensor[2]]])
    return Matrix


def gradToRad(angle):
    return (angle / 180) * np.pi


def Rx(xAngle):
    return np.array([[1, 0, 0],
                     [0, np.cos(xAngle), np.sin(xAngle)],
                     [0, -np.sin(xAngle), np.cos(xAngle)]])


def Ry(yAngle):
    return np.array([[np.cos(yAngle), 0, np.sin(yAngle)],
                     [0, 1, 0],
                     [np.sin(yAngle), 0, np.cos(yAngle)]])


def Rz(zAngle):
    return np.array([[np.cos(zAngle), np.sin(zAngle), 0],
                     [-np.sin(zAngle), np.cos(zAngle), 0],
                     [0, 0, 1]])


def RotFromTwoVecs(vec1, vec2):
    """ Find the rotation matrix that aligns vec1 to vec2
    :param vec1: A 3d "source" vector
    :param vec2: A 3d "destination" vector
    :return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
    """

    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
    return rotation_matrix


def radiusFromPt(pts, sigma):
    pts = np.abs(pts)
    if pts[1] > 0:
        teta = np.arctan(pts[2] / pts[1])
    else:
        teta = 0
    r = sigma[1] * sigma[2] / ((np.sin(teta) * sigma[2]) ** 2 + (np.cos(teta) * sigma[1]) ** 2) ** .5
    return r


def vecAbs(vec):
    return np.sqrt(vec[0] ** 2 + vec[1] ** 2 + vec[2] ** 2)


def vecDir(vec):
    return vec / vecAbs(vec)


def posVec(vec):
    return (vec ** 2) ** .5


def eulersFromRPG(R):
    tol = sys.float_info.epsilon * 10

    if abs(R.item(0, 0)) < tol and abs(R.item(1, 0)) < tol:
        eul1 = 0
        eul2 = m.atan2(-R.item(2, 0), R.item(0, 0))
        eul3 = m.atan2(-R.item(1, 2), R.item(1, 1))
    else:
        eul1 = m.atan2(R.item(1, 0), R.item(0, 0))
        sp = m.sin(eul1)
        cp = m.cos(eul1)
        eul2 = m.atan2(-R.item(2, 0), cp * R.item(0, 0) + sp * R.item(1, 0))
        eul3 = m.atan2(sp * R.item(0, 2) - cp * R.item(1, 2), cp * R.item(1, 1) - sp * R.item(0, 1))

    """
    print("z - phi =", eul1)
    print("y - theta  =", eul2)
    print("x - psi =", eul3)

    print("checkRPGRepo")
    print(R)

    print(np.dot(np.dot(Rz(eul1),Ry(eul2)),Rx(eul3)))
    """
    return eul1, eul2, eul3


def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = vecDir(v1)
    v2_u = vecDir(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))


def randomUnitVec():
    phi = np.random.uniform(0, np.pi * 2)
    costheta = np.random.uniform(-1, 1)

    theta = np.arccos(costheta)
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)
    return np.array([x, y, z])


def randomOrthMat():
    num_dim = 3
    x = special_ortho_group.rvs(num_dim)
    return x


def ellipsoidVol(sig):
    return 4 / 3 * np.pi * sig[0] * sig[1] * sig[2]


def calcAnisoMatrix(R):
    RPG = R  # np.linalg.eigh(R)[1]
    Rii = sum([RPG[0][0], RPG[1][1], RPG[2][2]])
    aniso = R / Rii - np.identity(3) / 3
    return aniso


def calcAnisoEigs(aniso):
    # eigenwerte für anisotropie muss hiermit berechnet werden!
    # wieso nochmal? es tauchen negative eigenwerte auf, bei berechnung mit numpy
    # nicht-symmetrische matrix (anisotrop)
    eigenVal = np.linalg.eig(aniso)
    eigens = list(eigenVal[0])
    eigVec = list(eigenVal[1])
    if list(eigens) == [0, 0, 0]:
        return np.zeros(3), None
    maxEigenIdx = eigens.index(max(eigens))
    minEigenIdx = eigens.index(min(eigens))
    middle = [0, 1, 2]
    middle.remove(maxEigenIdx)
    middle.remove(minEigenIdx)
    middle = middle[0]

    gamma_1 = eigens[maxEigenIdx]
    gamma_2 = eigens[middle]
    gamma_3 = eigens[minEigenIdx]

    eigenVal = np.array([gamma_1, gamma_2, gamma_3])
    return eigenVal, eigVec


def C_barycentric(R):
    aniso = calcAnisoMatrix(R)
    anisoEigs = calcAnisoEigs(aniso)[0]
    if list(anisoEigs) == [0, 0, 0]:
        return np.array([0, 0, 1])
    else:
        gamma_1 = anisoEigs[0]
        gamma_2 = anisoEigs[1]
        gamma_3 = anisoEigs[2]

    C1c = gamma_1 - gamma_2
    C2c = 2 * (gamma_2 - gamma_3)
    C3c = 3 * gamma_3 + 1
    CWeights = np.array([C1c, C2c, C3c])

    return CWeights
