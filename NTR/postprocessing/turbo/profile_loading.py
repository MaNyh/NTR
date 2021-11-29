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
        xl = bounds[1] - bounds[0]
        yl = bounds[3] - bounds[2]
        boxes.append(xl * yl)
    profile_id = boxes.index(min(boxes))
    profilepoints_ids = [idx for idx, i in enumerate(split["RegionId"]) if i == profile_id]
    profilepoints = split.extract_cells(profilepoints_ids)

    points, psPoly, ssPoly, ind_vk, ind_hk, midsPoly, metal_angle_vk, metal_angle_hk, camber_angle = extract_geo_paras(
        profilepoints.points, alpha)

    psVals = psPoly.sample(volmesh)
    ssVals = ssPoly.sample(volmesh)
    return psVals, ssVals, points, ind_vk, ind_hk, camber_angle


def calc_loading_volmesh(volmesh, alpha, verbose=False):
    psVals, ssVals, sortedPoints, ind_vk, ind_hk, camber_angle = extract_profile_from_volmesh(alpha, volmesh)

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

    ssVals.rotate_z(-camber_angle)
    psVals.rotate_z(-camber_angle)

    psVals = psVals.cell_data_to_point_data()
    ssVals = ssVals.cell_data_to_point_data()

    ps_xc = np.zeros(psVals.number_of_points)
    ps_cp = np.zeros(psVals.number_of_points)

    for idx, pt in enumerate(psVals.points):
        ps_xc[idx] = pt[0] / camber.length
        ps_cp[idx] = calc_inflow_cp(psVals.point_arrays["p"][idx], p2, p1)

    ss_xc = np.zeros(ssVals.number_of_points)
    ss_cp = np.zeros(ssVals.number_of_points)

    for idx, pt in enumerate(ssVals.points):
        ss_xc[idx] = pt[0] / camber.length
        ss_cp[idx] = calc_inflow_cp(ssVals.point_arrays["p"][idx], p2, p1)

    ssVals["xc"] = ss_xc
    psVals["xc"] = ps_xc

    if verbose:
        plt.figure()
        plt.plot(ss_xc, ss_cp)
        plt.plot(ps_xc, ps_cp)
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
    # todo: why are the first values here shitty?
    plt.figure()
    plt.plot(ss_xc, ss_cp, label="ss from experiment")
    plt.plot(ps_xc, ps_cp, label="ps from experiment")
    plt.plot(ss_xc_numerical[10:], ss_xc_cp_numerical[10:], label="ss from numeric")
    plt.plot(ps_xc_numerical[:-10], ps_xc_cp_numerical[:-10], label="ss from numeric")
    plt.show()
