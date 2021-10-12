# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 23:33:53 2020

@author: malte
"""

import os
import pyvista as pv
import vtk
import numpy as np

from NTR.utils.mathfunctions import vecAbs, unitVec, vecProjection, vecAngle
from NTR.utils.geom_functions.pyvista_utils import load_mesh
from NTR.utils.filehandling import yaml_dict_read



#################################################################################
#                                                                               #
#                                                                               #
#                          VectorHelper                                         #
#                                                                               #
#                                                                               #
#################################################################################


def cellDirections(cellUMean, wallNorm):
    x = unitVec(cellUMean)  # mainDirection
    z = unitVec(wallNorm)  # tangential
    y = unitVec(np.cross(x, z))  # spanwise
    return np.array([x, y, z])


#################################################################################
#                                                                               #
#                                                                               #
#                          VTK / VTKI Helper                                    #
#                                                                               #
#                                                                               #
#################################################################################


def constructWallMesh(meshList):
    wallMesh = pv.UnstructuredGrid()

    for mesh in meshList:
        m = load_mesh(mesh)
        wallMesh = wallMesh.merge(m)

    return wallMesh


def readDataSet(grid, dataName):
    grid.set_active_scalars(dataName)
    data = grid.active_scalars
    return data


def saveDataSets(grid, dataSets):
    for dataName, values in dataSets.items():
        grid[dataName] = values


def saveSolution(grid):
    # ÄNDERN in sinnvollen namen
    saveFileName = "gridSpacing.vtk"
    grid.save(saveFileName)


#################################################################################
#                                                                               #
#                                                                               #
#                           CalcFunctions                                       #
#                                                                               #
#                                                                               #
#################################################################################


def closestWallNormal(point,surfaceMesh):

    locator = vtk.vtkPointLocator()
    locator.SetDataSet(surfaceMesh)
    locator.BuildLocator()
    index = locator.FindClosestPoint(point)
    wallLocation = np.array(surfaceMesh.GetPoint(index))
    cellLocation = np.array(point)
    vector = wallLocation - cellLocation

    return vector


def calcWallNormalVectors(labelChunk,surfaceMesh):
    vectors = []

    for cellIdx in labelChunk:
        center = processData["cellCenters"][cellIdx]
        wallNormal = closestWallNormal(center, surfaceMesh)
        vectors.append(wallNormal)

    vectors = np.array(vectors)
    return vectors


def cellSpans(labelChunk, solutionMesh, processData):

    spans = []

    for cellIdx in labelChunk:

        x_span = 0
        x_weight = 0

        y_span = 0
        y_weight = 0

        z_span = 0
        z_weight = 0

        wallNormal = processData["wallNormal"][cellIdx]
        uMean = processData["UMean"][cellIdx]

        cellDirs = cellDirections(uMean, wallNormal)
        xx = cellDirs[0]
        yy = cellDirs[1]
        zz = cellDirs[2]

        egdeVectors = []

        CELL = solutionMesh.GetCell(cellIdx)
        edgeNumbers = CELL.GetNumberOfEdges()
        for edgeIdx in range(edgeNumbers):
            EDGE = CELL.GetEdge(edgeIdx)
            EDGE = EDGE.GetPoints()
            edgeVec = np.array(EDGE.GetPoint(0)) - np.array(EDGE.GetPoint(1))
            egdeVectors.append(edgeVec)

        for vec in egdeVectors:
            x_span += vecAbs(vecProjection(xx, vec))
            x_weight += abs(np.cos(vecAngle(xx, vec)))

            y_span += vecAbs(vecProjection(yy, vec))
            y_weight += abs(np.cos(vecAngle(yy, vec)))

            z_span += vecAbs(vecProjection(zz, vec))
            z_weight += abs(np.cos(vecAngle(zz, vec)))

        spans.append([x_span / x_weight, y_span / y_weight, z_span / z_weight])
        # pBarUpdate()
    return spans



def getWalluTaus(labelChunk, solutionMesh, mu_0, processData):

    uTaus = []
    for cellIdx in labelChunk:
        cellNormal = processData["wallNormal"][cellIdx]
        cellCenter = processData["cellCenters"][cellIdx]

        pointOnWall = cellCenter + cellNormal

        PT = pv.PolyData(pointOnWall)
        PT = PT.sample(solutionMesh)
        PT.set_active_scalars("gradient")

        gradUWall = PT.active_scalars[0]
        gradUWall = np.array([[gradUWall[0], gradUWall[1], gradUWall[2]], [gradUWall[3], gradUWall[4], gradUWall[5]],
                              [gradUWall[6], gradUWall[7], gradUWall[8]]])

        yDirection = unitVec(cellNormal)
        gradUyW = vecAbs(np.dot(yDirection, gradUWall))
        PT.set_active_scalars("rhoMean")
        rhoW = PT.active_scalars[0]
        tauW = gradUyW * mu_0 * rhoW
        u_tau = (tauW / rhoW) ** 0.5
        uTaus.append(u_tau)

    return uTaus


def gridSpacing(mu_0, processData):

    xSpans = processData["xSpan"]
    ySpans = processData["ySpan"]
    zSpans = processData["zSpan"]

    uTaus = processData["uTaus"]

    Deltax = xSpans * uTaus / mu_0
    Deltay = ySpans * uTaus / mu_0
    Deltaz = zSpans * uTaus / mu_0

    return [Deltax, Deltay, Deltaz]


#################################################################################
#                                                                               #
#                                                                               #
#                    construct DimlessGridSpacing                               #
#                                                                               #
#                                                                               #
#################################################################################

def calc(settings_yml):

    settings = yaml_dict_read(settings_yml)
    case_path = os.path.abspath(os.path.dirname(settings_yml))

    mu_0 = 2e-5

    processData = {}


    print("reading solutionMesh...")
    solutionMesh = load_mesh(solutionVTK)
    print("constructing surfacemesh from wall meshes ...")
    surfaceMesh = constructWallMesh(WallSurfacesVTKs)

    print("preparing processData from meshes")
    solutionMesh = solutionMesh.compute_derivative(scalars="UMean")
    processData["UMean"] = readDataSet(solutionMesh, "UMean")
    processData["cellCenters"] = solutionMesh.cell_centers().points

    cellIds = [i for i in range(solutionMesh.GetNumberOfCells())]

    print("calculating wall-normal vectors...")
    processData["wallNormal"] = calcWallNormalVectors(cellIds)

    print("calculating cell spans from WallNormals and CellEdges...")
    spanS = cellSpans(cellIds, surfaceMesh)
    processData["xSpan"] = np.array([i[0] for i in spanS])  # calculate cell span in flow direction
    processData["ySpan"] = np.array([i[1] for i in spanS])  # calculate cell span in wall normal direction
    processData["zSpan"] = np.array([i[2] for i in spanS])  # calculate cell span in span direction

    print("calculating wall-shear and friction-velocity")
    uTaus = getWalluTaus(cellIds, solutionMesh, mu_0, processData)
    processData["uTaus"] = uTaus

    print("calculating grid spacing")
    gridSpacings = gridSpacing(mu_0, processData)

    processData["DeltaXPlus"] = gridSpacings[0]
    print("min Dx+ : %.2f" % (min(processData["DeltaXPlus"])))
    print("max Dx+ : %.2f" % (max(processData["DeltaXPlus"])))

    processData["DeltaYPlus"] = gridSpacings[1]
    print("min Dy+ : %.2f" % (min(processData["DeltaYPlus"])))
    print("max Dy+ : %.2f" % (max(processData["DeltaYPlus"])))

    processData["DeltaZPlus"] = gridSpacings[2]
    print("min Dz+ : %.2f" % (min(processData["DeltaZPlus"])))
    print("max Dz+ : %.2f" % (max(processData["DeltaZPlus"])))

    print("calculation done, saving data to grid")
    solutionMesh = saveDataSets(solutionMesh, processData)
    saveSolution(solutionMesh)

