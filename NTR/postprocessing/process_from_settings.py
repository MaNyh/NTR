import os

from NTR.database.case_dirstructure import casedirs
from NTR.postprocessing.turbo.createProfileData import createProfileData
from NTR.utils.filehandling import yaml_dict_read, read_pickle
from NTR.utils.mathfunctions import vecAbs
from NTR.utils.mesh_handling.pyvista_utils import load_mesh
from NTR.utils.mesh_handling.pyvista_utils import mesh_scalar_gradients


def createProfileData_fromSettings(volmesh, settings):
    case_settings = yaml_dict_read(settings)
    case_path = os.path.dirname(settings)
    geo_path = os.path.join(case_path, casedirs["data"], "geometry.pkl")

    mplane_in = case_settings["post_settings"]["post_func_data"]["measureplane_slices"]["x_pos_1"]
    mplane_out = case_settings["post_settings"]["post_func_data"]["measureplane_slices"]["x_pos_2"]
    geomdat = read_pickle(geo_path)
    midspan_z = geomdat["span_z"] / 2
    alpha = case_settings["geometry"]["alpha"]
    # todo: insert a test-function for the definition of fluid-coeffs. incompressible sim? are values welldefined?
    kappa = float(case_settings["case_settings"]["fluid"]["kappa"])
    R_L = float(case_settings["case_settings"]["fluid"]["R_L"])
    p_k = float(case_settings["case_settings"]["fluid"]["p_k"])
    As = float(case_settings["case_settings"]["fluid"]["As"])
    cp = float(case_settings["case_settings"]["fluid"]["cp"])
    Ts = float(case_settings["case_settings"]["fluid"]["Ts"])
    l = vecAbs(
        geomdat["sortedPoly"][geomdat["hk_vk_idx"]["ind_vk"]] - geomdat["sortedPoly"][geomdat["hk_vk_idx"]["ind_hk"]])

    outputpath = os.path.join(case_path, casedirs["data"])
    createProfileData(volmesh, midspan_z, alpha, mplane_in, mplane_out, outputpath, kappa, R_L, p_k, As, l, cp, Ts)
