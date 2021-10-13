# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 15:53:51 2020

@author: malte
"""
import pyvista as pv
import numpy as np
import os
import matplotlib.pyplot as plt

#ToDo: this must be an argument of a function
meshPath = os.path.join("d:", "simulationsergebnisse", "tmpSolutions", "VTK", "1_Coarsest_462000.vtk")


def readVTKUnstructuredGrid(gridName):
    print("reading %s", gridName)
    meshData = pv.UnstructuredGrid(gridName)
    return meshData


def dataSetNames(grid):
    names = grid.array_names
    return names


def readDataSet(grid, name):
    grid.set_active_scalars(name)
    scalars = grid.active_scalars
    processData[name] = scalars


def readAllDataSets(grid):
    names = dataSetNames(grid)
    for name in names:
        readDataSet(grid, name)


def getSliceFromVTKGrid(grid, normalVec):
    vtkiSlice = grid.slice(normal=normalVec)
    return vtkiSlice


def calcRField(processData):
    resolvedStresses = processData["UPrime2Mean"]
    modelledStresses = processData["turbulenceProperties:RMean"]
    stressTensorField = resolvedStresses + modelledStresses

    Rlist = []
    for item in stressTensorField:
        uDash_uu = item[0]
        uDash_vv = item[1]
        uDash_ww = item[2]
        uDash_uv = item[3]
        uDash_vw = item[4]
        uDash_uw = item[5]
        R = np.array([[uDash_uu, uDash_uv, uDash_uw], [uDash_uv, uDash_vv, uDash_vw], [uDash_uw, uDash_vw, uDash_ww]])
        Rlist.append(R)

    return Rlist


def calckField(processData, RField):
    kList = []
    for R in RField:
        k = 0.5 * (R[0][0] + R[1][1] + R[2][2])
        kList.append(k)

    kField = np.array(kList)
    return kField


def calcAnisoMatrix(processData, RField, kField):
    anisoList = []
    for idx in range(len(RField)):
        aniso = RField[idx] / (2 * kField[idx]) - np.identity(3) / 3
        anisoList.append(aniso)
    anisoField = np.array(anisoList)

    return anisoField


def eigenValAniso(AnisoField):
    eigenValField = []
    # diagonalField = []
    for aniso in AnisoField:
        eigenVal = np.linalg.eigh(
            aniso)  # eigh für symmetrische matritzen --> aaber dann ist die diagonale nicht mehr sortiert als einheitsmatrix
        eigens = list(eigenVal[0])
        maxEigenIdx = eigens.index(max(eigens))
        minEigenIdx = eigens.index(min(eigens))
        middle = [0, 1, 2]
        middle.remove(maxEigenIdx)
        middle.remove(minEigenIdx)
        middle = middle[0]

        gamma_1 = eigens[maxEigenIdx]
        gamma_2 = eigens[middle]
        gamma_3 = eigens[minEigenIdx]

        eigenValField.append(np.array([gamma_1, gamma_2, gamma_3]))

    return eigenValField


def II_III(eigenValField):
    IIField = []
    IIIField = []
    for eigenVals in eigenValField:
        gamma_1 = eigenVals[0]
        gamma_2 = eigenVals[1]
        II = gamma_1 ** 2 + gamma_1 * gamma_2 + gamma_2 ** 2
        III = -gamma_1 * gamma_2 * (gamma_1 + gamma_2)
        IIField.append(II)
        IIIField.append(III)
    return IIField, IIIField


def lumleyTriangle(III, II):
    # fig=plt.figure()
    plt.vlines(0, min(II), max(II), linestyle='--')
    plt.scatter(III, II, marker='.', label="solution")

    ############################################
    #                                          #
    #                1C, 2C, 3C                #
    #                                          #
    ############################################

    eigField1C = [[2 / 3, -1 / 3, -1 / 3]]
    C1_II, C1_III = II_III(eigField1C)

    eigField2C = [[1 / 6, 1 / 6, -1 / 3]]
    C2_II, C2_III = II_III(eigField2C)

    eigField3C = [[0, 0, 0]]
    C3_II, C3_III = II_III(eigField3C)

    plt.scatter(C1_III, C1_II, label='XC1')
    plt.scatter(C2_III, C2_II, label='XC2')
    plt.scatter(C3_III, C3_II, label='XC3')

    ############################################
    #                                          #
    #                1C  to  3C                #
    #                                          #
    ############################################

    res = 200
    C1_to_C3_gamma_1 = np.linspace(2 / 3, 0, res, endpoint=False)
    C1_to_C3_gamma_2 = np.linspace(-1 / 3, 0, res, endpoint=False)
    eigFieldC1toC2 = [[C1_to_C3_gamma_1[i], C1_to_C3_gamma_2[i], 0] for i in range(res)]
    C1to2_II, C1to2_III = II_III(eigFieldC1toC2)
    plt.plot(C1to2_III, C1to2_II, color='k')

    ############################################
    #                                          #
    #                1C  to  2C                #
    #                                          #
    ############################################

    C1_to_C2_gamma_1 = np.linspace(2 / 3, 1 / 6, res, endpoint=False)
    C1_to_C2_gamma_2 = np.linspace(-1 / 3, 1 / 6, res, endpoint=False)
    eigFieldC1toC2 = [[C1_to_C2_gamma_1[i], C1_to_C2_gamma_2[i], 0] for i in range(res)]
    C1to2_II, C1to2_III = II_III(eigFieldC1toC2)
    plt.plot(C1to2_III, C1to2_II, color='r')

    ############################################
    #                                          #
    #                3C  to  2C                #
    #                                          #
    ############################################

    C3_to_C2_gamma_1 = np.linspace(0, 1 / 6, res, endpoint=False)
    C3_to_C2_gamma_2 = np.linspace(0, 1 / 6, res, endpoint=False)
    eigFieldC1toC2 = [[C3_to_C2_gamma_1[i], C3_to_C2_gamma_2[i], 0] for i in range(res)]
    C1to2_II, C1to2_III = II_III(eigFieldC1toC2)
    plt.plot(C1to2_III, C1to2_II, color='b')

    plt.title("Lumley triangle")
    plt.legend()
    plt.show()
    return C1_II, C3_II


def eigenValueMap(gamma_2, gamma_1):
    plt.vlines(0, 0, 0.7, linestyle='--')
    plt.scatter(gamma_2, gamma_1, marker='.')
    plt.title("Eigenvalue map")

    ############################################
    #                                          #
    #                1C, 2C, 3C                #
    #                                          #
    ############################################

    eigField1C = [2 / 3, -1 / 3, -1 / 3]
    eigField2C = [1 / 6, 1 / 6, -1 / 3]
    eigField3C = [0, 0, 0]

    plt.scatter(eigField1C[1], eigField1C[0], label='XC1')
    plt.scatter(eigField2C[1], eigField2C[0], label='XC2')
    plt.scatter(eigField3C[1], eigField3C[0], label='XC3')

    ############################################
    #                                          #
    #                1C  to  3C                #
    #                                          #
    ############################################
    res = 200
    C1_to_C3_gamma_1 = np.linspace(2 / 3, 0, res, endpoint=False)
    C1_to_C3_gamma_2 = np.linspace(-1 / 3, 0, res, endpoint=False)
    plt.plot(C1_to_C3_gamma_2, C1_to_C3_gamma_1, color='k')

    ############################################
    #                                          #
    #                1C  to  2C                #
    #                                          #
    ############################################
    C1_to_C2_gamma_1 = np.linspace(2 / 3, 1 / 6, res, endpoint=False)
    C1_to_C2_gamma_2 = np.linspace(-1 / 3, 1 / 6, res, endpoint=False)
    plt.plot(C1_to_C2_gamma_2, C1_to_C2_gamma_1, color='r')

    ############################################
    #                                          #
    #                3C  to  1C                #
    #                                          #
    ############################################

    C3_to_C2_gamma_1 = np.linspace(0, 1 / 6, res, endpoint=False)
    C3_to_C2_gamma_2 = np.linspace(0, 1 / 6, res, endpoint=False)
    plt.plot(C3_to_C2_gamma_2, C3_to_C2_gamma_1, color='b')

    plt.legend()
    plt.show()


def turbTriangle(Eta, Xi):
    plt.scatter(Xi, Eta)
    plt.vlines(0, 0, 0.3, linestyle='--')
    ############################################
    #                                          #
    #                1C  to  3C                #
    #                                          #
    ############################################

    res = 200
    C1_to_C3_gamma_1 = np.linspace(2 / 3, 0, res, endpoint=False)
    C1_to_C3_gamma_2 = np.linspace(-1 / 3, 0, res, endpoint=False)
    eigFieldC1toC3 = [[C1_to_C3_gamma_1[i], C1_to_C3_gamma_2[i], 0] for i in range(res)]
    C1to3_II, C1to3_III = II_III(eigFieldC1toC3)
    Eta_1to3 = np.sqrt(np.array(C1to3_II) / 3)
    Xi_1to3 = np.cbrt(np.array(C1to3_III) / 2)
    plt.plot(Xi_1to3, Eta_1to3, color='k')

    ############################################
    #                                          #
    #                1C  to  2C                #
    #                                          #
    ############################################

    C1_to_C2_gamma_1 = np.linspace(2 / 3, 1 / 6, res, endpoint=False)
    C1_to_C2_gamma_2 = np.linspace(-1 / 3, 1 / 6, res, endpoint=False)
    eigFieldC1toC2 = [[C1_to_C2_gamma_1[i], C1_to_C2_gamma_2[i], 0] for i in range(res)]
    C1to2_II, C1to2_III = II_III(eigFieldC1toC2)
    Eta_1to2 = np.sqrt(np.array(C1to2_II) / 3)
    Xi_1to2 = np.cbrt(np.array(C1to2_III) / 2)
    plt.plot(Xi_1to2, Eta_1to2, color='r')

    ############################################
    #                                          #
    #                3C  to  2C                #
    #                                          #
    ############################################

    C3_to_C2_gamma_1 = np.linspace(0, 1 / 6, res, endpoint=False)
    C3_to_C2_gamma_2 = np.linspace(0, 1 / 6, res, endpoint=False)
    eigFieldC3toC2 = [[C3_to_C2_gamma_1[i], C3_to_C2_gamma_2[i], 0] for i in range(res)]
    C3to2_II, C3to2_III = II_III(eigFieldC3toC2)

    Eta_3to2 = np.sqrt(np.array(C3to2_II) / 3)
    Xi_3to2 = np.cbrt(np.array(C3to2_III) / 2)
    plt.plot(Xi_3to2, Eta_3to2, color='b')

    plt.title("turbulence triangle")
    plt.legend()
    plt.show()


"""
def barycentricMap(eigenValueField):

    XC1 = np.array([1,0])
    XC2 = np.array([0,0])
    XC3 = np.array([1/2,np.cbrt(3)/2])

    plotx= []
    ploty= []

    for eigenVals in eigenValueField:
        gamma_1 = eigenVals[0]
        gamma_2 = eigenVals[1]
        gamma_3 = eigenVals[2]

        C1c = gamma_1-gamma_2
        C2c = 2*(gamma_2-gamma_1)
        C3c = 3*(gamma_3)+1

    plt.legend()
    plt.title("barycentric map")
    plt.show()
    return 0
"""

if __name__ == "__main__":
    processData = {}
    calcedData = {}

    meshObj = readVTKUnstructuredGrid(meshPath)
    sliceObj = getSliceFromVTKGrid(meshObj, [0, 0, 1])
    sliceObj.set_active_scalars("UPrime2Mean")
    sliceObj.plot(show_edges=True, cpos=(0, 0, 1))

    # Slice or complete Solution?
    readAllDataSets(meshObj)

    RField = calcRField(processData)
    calcedData["RField"] = RField
    kField = calckField(processData, RField)
    calcedData["kField"] = kField
    AnisoField = calcAnisoMatrix(processData, RField, kField)
    calcedData["AnisoField"] = AnisoField
    EigenValues = eigenValAniso(AnisoField)
    calcedData["EigenValues"] = EigenValues
    InvII, InvIII = II_III(EigenValues)
    calcedData["InvII"] = InvII
    calcedData["InvIII"] = InvIII

    eta = np.sqrt(np.array(InvII) / 3)
    xi = np.cbrt(np.array(InvIII) / 2)

    lumleyTriangle(InvIII, InvII)
    eigenValueMap([i[1] for i in EigenValues], [i[0] for i in EigenValues])
    turbTriangle(eta, xi)

