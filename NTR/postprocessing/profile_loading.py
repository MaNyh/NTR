import os
import tempfile
import pyvista as pv
import numpy as np
from matplotlib import pyplot as plt

from NTR.utils.externals.paraview.cgns_to_vtk import convert_cgns_to_vtk
from NTR.utils.geom_functions.pyvista_utils import load_mesh
from NTR.utils.geom_functions.distance import closest_node_index
from NTR.preprocessing.create_geom import extract_geo_paras
from NTR.utils.fluid_functions.aeroFunctions import calc_inflow_cp
from NTR.utils.filehandling import yaml_dict_read
from NTR.database.case_dirstructure import casedirs

settings_yaml_file = "D:\\NTR\\examples\\CascadeCase_gwk_rans_trace\\case_settings.yml"


def extract_profile_from_volmesh(settings_yml, volmesh):
    settings = yaml_dict_read(settings_yml)
    alpha = settings["geometry"]["alpha"]
    bounds = volmesh.bounds
    midspan_z = (bounds[-1] - bounds[-2]) / 2
    z_slice = volmesh.slice(normal="z", origin=(0, 0, midspan_z))
    edges = z_slice.extract_feature_edges()
    split = edges.connectivity()
    profilepoints_ids = [idx for idx, i in enumerate(split["RegionId"]) if i == 0]
    profilepoints = split.extract_cells(profilepoints_ids)

    points, psPoly, ssPoly, ind_vk, ind_hk, midsPoly, metal_angle_vk, metal_angle_hk, camber_angle = extract_geo_paras(
        profilepoints.points, alpha)

    pspts = []
    for pt in psPoly.points:
        closest = closest_node_index(pt, profilepoints.points)
        pspts.append(closest)
    sspts = []
    for pt in ssPoly.points:
        closest = closest_node_index(pt, profilepoints.points)
        sspts.append(closest)

    psVals = profilepoints.extract_points(pspts)
    ssVals = profilepoints.extract_points(sspts)
    return psVals, ssVals, points, ind_vk, ind_hk, camber_angle

def paraview_convert_cgns_to_vtk(settings_yml):
    settings = yaml_dict_read(settings_yml)
    case_path = os.path.dirname(settings_yml)

    assert "post_settings" in settings.keys(), "post_settings missing"
    assert "paraview_convert_cgns_to_vtk" in settings["post_settings"].keys(), "paraview_convert_cgns_to_vtk settings missing"
    assert "cgns" in settings["post_settings"]["paraview_convert_cgns_to_vtk"].keys(), "no cgns to convert defined in settings"
    path_to_volmesh = os.path.join(case_path, settings["post_settings"]["paraview_convert_cgns_to_vtk"]["cgns"])

    assert os.path.isfile(path_to_volmesh), "file " + path_to_volmesh + " does not exist"

    target_name = settings["post_settings"]["paraview_convert_cgns_to_vtk"]["vtk_name"]
    convert_cgns_to_vtk(path_to_volmesh, os.path.join(case_path, casedirs["solution"], target_name))
    return 0

def calc_loading_volmesh(settings_yml):
    settings = yaml_dict_read(settings_yml)
    case_path = os.path.dirname(settings_yml)

    path_to_volmesh = os.path.join(case_path, settings["post_settings"]["paraview_convert_cgns_to_vtk"]["cgns"])

    assert os.path.isfile(path_to_volmesh), "file " + path_to_volmesh + " does not exist"
    temp_dir = tempfile.TemporaryDirectory()

    mesh_ext = os.path.splitext(path_to_volmesh)[-1]
    volmesh = None

    if mesh_ext == ".cgns":
        tmpmesh = "solution.vtk"
        convert_cgns_to_vtk(path_to_volmesh, os.path.join(case_path, casedirs["solution"], tmpmesh))
        volmesh = load_mesh(os.path.join(temp_dir.name, os.path.join(case_path, casedirs["solution"], tmpmesh)))
        temp_dir.cleanup()
    elif mesh_ext == ".vtk":
        volmesh(path_to_volmesh)

    psVals, ssVals, points, ind_vk, ind_hk, camber_angle = extract_profile_from_volmesh(settings_yml, volmesh)

    return psVals, ssVals, points, ind_vk, ind_hk, camber_angle


"""
ssVals, psVals, sortedPoints, ind_vk, ind_hk, camber_angle = calc_loading_volmesh(settings_yaml_file)

camber = pv.Line((0, 0, 0), -(sortedPoints[ind_vk] - sortedPoints[ind_hk]))
xLine = pv.Line((-1, 0, 0), (1, 0, 0))
yLine = pv.Line((0, -1, 0), (0, 1, 0))


shift = sortedPoints[ind_vk]
shift -= psVals.points[0][-1]

ssVals.points -= shift
psVals.points -= shift

ssVals.rotate_z(-camber_angle + 90)
psVals.rotate_z(-camber_angle + 90)

psVals = psVals.cell_data_to_point_data()
ssVals = ssVals.cell_data_to_point_data()

ps_xc = np.zeros(ssVals.number_of_points)
ps_cp = np.zeros(ssVals.number_of_points)

for idx, pt in enumerate(psVals.points):
    ps_xc[idx] = pt[0] / camber.length
    ps_cp[idx] = calc_inflow_cp(psVals.point_arrays["Pressure"][idx],101315,98000)

psVals["xc"] = ps_xc

ss_xc = np.zeros(ssVals.number_of_points)
ss_cp = np.zeros(psVals.number_of_points)

for idx, pt in enumerate(ssVals.points):
    ss_xc[idx] = pt[0] / camber.length

ssVals["xc"] = ss_xc

p = pv.Plotter()
p.add_mesh(camber)
p.add_mesh(ssVals)
p.add_mesh(psVals)
p.add_mesh(xLine, color="yellow")
p.add_mesh(yLine)
p.show()

psVals_monotone_xc = ps_xc#np.maximum.accumulate(psVals["xc"])
ssVals_monotone_xc = ss_xc#np.maximum.accumulate(ssVals["xc"])

fig = plt.figure()
plt.plot(ssVals_monotone_xc, ssVals.point_arrays["Pressure"])
plt.plot(psVals_monotone_xc, psVals.point_arrays["Pressure"])
plt.show()
"""
paraview_convert_cgns_to_vtk(settings_yaml_file)

