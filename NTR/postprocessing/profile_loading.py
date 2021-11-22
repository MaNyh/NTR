import os

import pyvista as pv
import numpy as np
from matplotlib import pyplot as plt

from NTR.database.datasets.measurement_data_2016111_gwk_compressor_gilge.interpret_raw import read_gilgegwk
from NTR.postprocessing.sim_values import getXSliceVals
from NTR.utils.filehandling import yaml_dict_read
from NTR.utils.fluid_functions.aeroFunctions import calc_inflow_cp
from NTR.preprocessing.create_geom import extract_geo_paras
from NTR.utils.mesh_handling.pyvista_utils import load_mesh
from NTR.utils.mathfunctions import absvec_array

def extract_profile_from_volmesh(alpha, volmesh):
    bounds = volmesh.bounds
    midspan_z = (bounds[-1] - bounds[-2]) / 2
    z_slice = volmesh.slice(normal="z", origin=(0, 0, midspan_z))

    edges = z_slice.extract_feature_edges()
    split = edges.connectivity()

    regionIds = list(dict.fromkeys(split["RegionId"]))
    boxes = []
    for rId in regionIds:
        regionPoints_ids = [idx for idx, i in enumerate(split["RegionId"]) if i == rId]
        regionPoints = split.extract_cells(regionPoints_ids)
        bounds = regionPoints.bounds
        xl = bounds[1]-bounds[0]
        yl = bounds[3]-bounds[2]
        boxes.append(xl*yl)
    profile_id = boxes.index(min(boxes))
    profilepoints_ids = [idx for idx, i in enumerate(split["RegionId"]) if i == profile_id]
    profilepoints = split.extract_cells(profilepoints_ids)

    points, psPoly, ssPoly, ind_vk, ind_hk, midsPoly, metal_angle_vk, metal_angle_hk, camber_angle = extract_geo_paras(
        profilepoints.points, alpha)

    psVals = psPoly.sample(volmesh)
    ssVals = ssPoly.sample(volmesh)
    return psVals, ssVals, points, ind_vk, ind_hk, camber_angle


def calc_loading_volmesh(volmesh,alpha, verbose=False):
    psVals, ssVals, sortedPoints, ind_vk, ind_hk, camber_angle = extract_profile_from_volmesh(alpha, volmesh)

    camber = pv.Line((0, 0, 0), -(sortedPoints[ind_vk] - sortedPoints[ind_hk]))
    #xLine = pv.Line((-1, 0, 0), (1, 0, 0))
    #yLine = pv.Line((0, -1, 0), (0, 1, 0))

    shift = sortedPoints[ind_vk]
    shift -= psVals.points[0][-1]

    ssVals.points -= shift
    psVals.points -= shift

    ssVals.rotate_z(-camber_angle + 90)
    psVals.rotate_z(-camber_angle + 90)

    psVals = psVals.cell_data_to_point_data()
    ssVals = ssVals.cell_data_to_point_data()

    ps_xc = np.zeros(psVals.number_of_points)
    ps_cp = np.zeros(psVals.number_of_points)

    for idx, pt in enumerate(psVals.points):
        ps_xc[idx] = pt[0] / camber.length
        ps_cp[idx] = calc_inflow_cp(psVals.point_arrays["Pressure"][idx], 101315, 98000)

    ss_xc = np.zeros(ssVals.number_of_points)
    ss_cp = np.zeros(ssVals.number_of_points)

    for idx, pt in enumerate(ssVals.points):
        ss_xc[idx] = pt[0] / camber.length
        ss_cp[idx] = calc_inflow_cp(ssVals.point_arrays["Pressure"][idx], 101315, 98000)


    ssVals["xc"] = ss_xc
    psVals["xc"] = ps_xc

    ssVals_monotone_xc = []
    ssVals_monotone_Pressure = []

    for idx in range(ssVals.number_of_points):
        if len(ssVals_monotone_xc) == 0 or ssVals["xc"][idx] > ssVals_monotone_xc[-1]:
            ssVals_monotone_xc.append(ssVals["xc"][idx])
            ssVals_monotone_Pressure.append(ssVals["Pressure"][idx])

    psVals_monotone_xc = []
    psVals_monotone_Pressure = []

    for idx in range(psVals.number_of_points):
        if len(psVals_monotone_xc) == 0 or psVals["xc"][idx] > psVals_monotone_xc[-1]:
            psVals_monotone_xc.append(psVals["xc"][idx])
            psVals_monotone_Pressure.append(psVals["Pressure"][idx])

    if verbose:
        plt.figure()
        plt.plot(ssVals_monotone_xc, ssVals_monotone_Pressure)
        plt.plot(psVals_monotone_xc, psVals_monotone_Pressure)
        plt.show()
    return psVals, ssVals


def compare_profileloading_numexp(settings_yml):

    settings = yaml_dict_read(settings_yml)
    case_path = os.path.dirname(settings_yml)
    path_to_volmesh = os.path.join(case_path, settings["post_settings"]["volmesh"])
    assert os.path.isfile(path_to_volmesh), "file " + path_to_volmesh + " does not exist"

    volmesh = load_mesh(path_to_volmesh)
    alpha = settings["geometry"]["alpha"]
    psVals, ssVals = calc_loading_volmesh(volmesh, alpha)

    inlet = getXSliceVals(volmesh, 0)
    inlet_sizes = inlet.compute_cell_sizes()

    if "Velocity" in inlet_sizes.array_names:
        U = inlet_sizes["Velocity"]
    elif "UMean" in inlet_sizes.array_names:
        U = inlet_sizes["UMean"]
    if "Density" in inlet_sizes.array_names:
        rho = inlet_sizes["Density"]
    elif "rhoMean" in inlet_sizes.array_names:
        rho = inlet_sizes["rhoMean"]
    if "Pressure" in inlet_sizes.array_names:
        p = inlet_sizes["Pressure"]
    elif "pMean" in inlet_sizes.array_names:
        p = inlet_sizes["pMean"]

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
