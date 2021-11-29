import os

from NTR.database.case_dirstructure import casedirs
from NTR.postprocessing.turbo.createProfileData import createProfileData
from NTR.postprocessing.turbo.profile_loading import calc_loading_volmesh
from NTR.utils.filehandling import yaml_dict_read, read_pickle
from NTR.utils.mathfunctions import vecAbs
from NTR.utils.fluid_functions.fluids import idealgas

def createProfileData_fromSettings(volmesh, settings):
    case_settings = yaml_dict_read(settings)
    case_path = os.path.dirname(settings)
    geo_path = os.path.join(case_path, casedirs["data"], "geometry.pkl")

    mplane_in = case_settings["post_settings"]["post_func_data"]["measureplane_slices"]["x_pos_1"]
    mplane_out = case_settings["post_settings"]["post_func_data"]["measureplane_slices"]["x_pos_2"]
    geomdat = read_pickle(geo_path)
    midspan_z = geomdat["span_z"] / 2
    alpha = case_settings["geometry"]["alpha"]

    fluid = idealgas(**case_settings["case_settings"]["fluid"])
    kappa = fluid.kappa
    Rs = fluid.Rs
    cp = fluid.cp

    #p_k = float(case_settings["post_settings"]["post_func_data"]["p_k"])
    As = float(case_settings["post_settings"]["post_func_data"]["As"])
    Ts = float(case_settings["post_settings"]["post_func_data"]["Ts"])

    l = vecAbs(
        geomdat["sortedPoly"][geomdat["hk_vk_idx"]["ind_vk"]] - geomdat["sortedPoly"][geomdat["hk_vk_idx"]["ind_hk"]])

    outputpath = os.path.join(case_path, casedirs["data"])
    return createProfileData(volmesh, midspan_z, alpha, mplane_in, mplane_out, outputpath, kappa, Rs, p_k, As, l, cp, Ts)


def computeProfileLoading_fromSettings(volmesh,settings):
    case_settings = yaml_dict_read(settings)
    alpha = case_settings["geometry"]["alpha"]
    return calc_loading_volmesh(volmesh, alpha)
