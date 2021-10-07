import os
import tempfile
import numpy as np

from NTR.utils.externals.paraview.cgns_to_vtk import convert_cgns_to_vtk
from NTR.utils.geom_functions.pyvista_utils import load_mesh
from NTR.preprocessing.create_geom import extract_geo_paras
from NTR.utils.filehandling import yaml_dict_read
from NTR.database.case_dirstructure import casedirs

settings_yaml_file = "D:\\NTR\\examples\\CascadeCase_gwk_rans_trace\\case_settings.yml"


def extract_profile_from_volmesh(settings_yml,volmesh):
    settings = yaml_dict_read(settings_yml)
    alpha = settings["geometry"]["alpha"]
    bounds = volmesh.bounds
    midspan_z = (bounds[-1] - bounds[-2])/2
    z_slice = volmesh.slice(normal="z",origin=(0,0,midspan_z))
    edges = z_slice.extract_feature_edges()
    split = edges.connectivity()
    profilepoints_ids = [idx for idx, i in enumerate(split["RegionId"]) if i==0]
    profilepoints = split.extract_cells(profilepoints_ids)

    points, psPoly, ssPoly, ind_vk, ind_hk, midsPoly, metal_angle_vk, metal_angle_hk, camber_angle = extract_geo_paras(profilepoints.points, alpha)

    for pt in psPoly.points:
        np.where(pt,profilepoints)
    return 0


def calc_loading_volmesh(settings_yml):
    settings = yaml_dict_read(settings_yml)
    case_path = os.path.dirname(settings_yml)

    path_to_volmesh = os.path.join(case_path,settings["solution"]["volmesh"])

    assert os.path.isfile(path_to_volmesh), "file " + path_to_volmesh + " does not exist"
    temp_dir = tempfile.TemporaryDirectory()

    mesh_ext = os.path.splitext(path_to_volmesh)[-1]
    volmesh = None

    if mesh_ext == ".cgns":
        tmpmesh = "solution.vtk"
        convert_cgns_to_vtk(path_to_volmesh, os.path.join(case_path,casedirs["solution"],tmpmesh))
        volmesh = load_mesh(os.path.join(temp_dir.name, os.path.join(case_path,casedirs["solution"],tmpmesh)))
        temp_dir.cleanup()
    elif mesh_ext == ".vtk":
        volmesh(path_to_volmesh)

    profilesolution = extract_profile_from_volmesh(settings_yml,volmesh)
    return profilesolution


profilesolution = calc_loading_volmesh(settings_yaml_file)
