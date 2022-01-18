import os
import pyvista as pv
import numpy as np
from matplotlib import pyplot as plt

from NTR.database.datasets.measurement_data_2016111_gwk_compressor_gilge.interpret_raw import read_gilgegwk
from NTR.postprocessing.generic.sim_values import getXSliceVals
from NTR.utils.filehandling import yaml_dict_read
from NTR.utils.fluid_functions.aeroFunctions import calc_inflow_cp
from NTR.preprocessing.create_geom import extract_geo_paras
from NTR.utils.mesh_handling.pyvista_utils import load_mesh
from NTR.utils.mathfunctions import absvec_array, absVec
from NTR.utils.mathfunctions import vecProjection
from NTR.utils.mesh_handling.pyvista_utils import lines_from_points
from NTR.preprocessing.create_geom import calcConcaveHull

def extract_profile_from_volmesh(alpha, volmesh):
    surf_idx = volmesh.surface_indices()
    surface_vol = volmesh.extract_cells(surf_idx)
    surface = surface_vol.extract_surface()
    bounds = surface.bounds

    midspan_z = (bounds[-1] - bounds[-2]) / 2
    z_slice = surface.slice(normal="z", origin=(0, 0, midspan_z))
    z_slice["pointIds"] = [i for i in range(z_slice.number_of_points)]

    z_slicepoints = z_slice.points

    allxx, allyy = z_slicepoints[::,0],z_slicepoints[::,1]
    outerxx, outeryy = calcConcaveHull(z_slicepoints[::,0],z_slicepoints[::,1],alpha)

    outerIds = []
    innerIds = []
    for idx in z_slice["pointIds"]:
        if allxx[idx] in list(outerxx) and allyy[idx] in list(outeryy):
            outerIds.append(idx)
        else:
            innerIds.append(idx)

    profilepoints = z_slice.extract_points(innerIds)


    points, psPoly, ssPoly, ind_vk, ind_hk, midsPoly, metal_angle_vk, metal_angle_hk, camber_angle = extract_geo_paras(
        profilepoints.points, alpha)
    psVals = lines_from_points(psPoly.points)
    ssVals = lines_from_points(ssPoly.points)
    psVals = psVals.sample(surface)
    ssVals = ssVals.sample(surface)
    ssBladeForce = 0
    psBladeForce = 0

    for i in range(psVals.number_of_cells):
        cell = psVals.extract_cells(i)
        pc = cell["p"][0]
        cellVec = cell.points[1]-cell.points[0]
        cellProject = vecProjection((1,0,0),cellVec)
        projectedLength = absVec(cellProject)
        yforce = pc*projectedLength
        psBladeForce += yforce

    for i in range(ssVals.number_of_cells):
        cell = ssVals.extract_cells(i)
        pc = cell["p"][0]
        cellVec = cell.points[1] - cell.points[0]
        cellProject = vecProjection((1, 0, 0), cellVec)
        projectedLength = absVec(cellProject)
        yforce = pc * projectedLength
        ssBladeForce += yforce

    bladeForce = psBladeForce - ssBladeForce

    return psVals, ssVals, points, ind_vk, ind_hk, camber_angle, bladeForce


def calc_loading_volmesh(volmesh, alpha, verbose=False):
    psVals, ssVals, sortedPoints, ind_vk, ind_hk, camber_angle, bladeForce_1d = extract_profile_from_volmesh(alpha, volmesh)

    bounds = volmesh.bounds
    x1 = bounds[0] + 1e-5 * bounds[1]
    x2 = bounds[1] + 1e-5 * bounds[0]

    inlet = volmesh.slice(normal="x", origin=(x1, 0, 0)).compute_cell_sizes().point_data_to_cell_data()
    outlet = volmesh.slice(normal="x", origin=(x2, 0, 0)).compute_cell_sizes().point_data_to_cell_data()
    p1 = sum(inlet["p"] * inlet["Area"]) / sum(inlet["Area"])
    p2 = sum(outlet["p"] * outlet["Area"]) / sum(outlet["Area"])

    camber = pv.Line((0, 0, 0), -(sortedPoints[ind_vk] - sortedPoints[ind_hk]))

    shift = sortedPoints[ind_vk]
    shift -= psVals.points[0][-1]

    ssVals.points -= shift
    psVals.points -= shift

    ssVals.rotate_z(-camber_angle,inplace=True)
    psVals.rotate_z(-camber_angle,inplace=True)

    psVals = psVals.cell_data_to_point_data()
    ssVals = ssVals.cell_data_to_point_data()

    ps_xc = np.zeros(psVals.number_of_points)
    ps_cp = np.zeros(psVals.number_of_points)

    for idx, pt in enumerate(psVals.points):
        ps_xc[idx] = pt[0] / camber.length
        ps_cp[idx] = calc_inflow_cp(psVals.point_data["p"][idx], p2, p1)

    ss_xc = np.zeros(ssVals.number_of_points)
    ss_cp = np.zeros(ssVals.number_of_points)

    for idx, pt in enumerate(ssVals.points):
        ss_xc[idx] = pt[0] / camber.length
        ss_cp[idx] = calc_inflow_cp(ssVals.point_data["p"][idx], p2, p1)

    ssVals["xc"] = ss_xc
    ssVals["cp"] = ss_cp
    psVals["xc"] = ps_xc
    psVals["cp"] = ps_cp

    if verbose:
        plt.figure()
        plt.plot(ss_xc, ss_cp)
        plt.plot(ps_xc, ps_cp)
        plt.show()
    return psVals, ssVals ,bladeForce_1d


def compare_profileloading_numexp(settings_yml):
    settings = yaml_dict_read(settings_yml)
    case_path = os.path.dirname(settings_yml)
    path_to_volmesh = os.path.join(case_path, settings["post_settings"]["volmesh"])
    assert os.path.isfile(path_to_volmesh), "file " + path_to_volmesh + " does not exist"

    volmesh = load_mesh(path_to_volmesh)
    alpha = settings["geometry"]["alpha"]
    psVals, ssVals, bladeForce = calc_loading_volmesh(volmesh, alpha)

    inlet = getXSliceVals(volmesh, 0)
    inlet_sizes = inlet.compute_cell_sizes()

    pdyn = absvec_array(U) ** 2 * rho / 2

    pressure = sum(p * inlet_sizes["Area"]) / sum(inlet_sizes["Area"])
    pressure_tot = sum((p + pdyn) * inlet_sizes["Area"]) / sum(inlet_sizes["Area"])

    ps_xc_numerical = psVals["xc"]
    ps_xc_pressure_numerical = psVals["Pressure"]
    ps_xc_cp_numerical = calc_inflow_cp(ps_xc_pressure_numerical, pressure_tot, pressure)
    ss_xc_numerical = ssVals["xc"]
    ss_xc_pressure_numerical = ssVals["Pressure"]
    ss_xc_cp_numerical = calc_inflow_cp(ss_xc_pressure_numerical, pressure_tot, pressure)
    ss_xc, ss_cp, ps_xc, ps_cp = read_gilgegwk(verbose=False)

    plt.figure()
    plt.plot(ss_xc, ss_cp, label="ss from experiment")
    plt.plot(ps_xc, ps_cp, label="ps from experiment")
    plt.plot(ss_xc_numerical, ss_xc_cp_numerical, label="ss from numeric")
    plt.plot(ps_xc_numerical, ps_xc_cp_numerical, label="ss from numeric")
    plt.show()
