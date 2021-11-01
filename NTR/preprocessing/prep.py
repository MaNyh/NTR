import pyvista as pv
import os
import numpy as np

from NTR.utils.filehandling import yaml_dict_read, write_pickle
from NTR.database.case_dirstructure import casedirs


def prep_channelflow_geometry(settings_yml):
    settings = yaml_dict_read(settings_yml)

    length = float(settings["geometry"]["channel_length"])
    width = float(settings["geometry"]["yper_width"])
    halfheight = float(settings["geometry"]["channel_halfheight"])

    yper_low = pv.Line((0, 0, 0), (length, 0, 0))
    yper_high = pv.Line((0, 0, width), (length, 0, width))
    span_z = halfheight

    x_mcl = np.arange(0, length, length / 100)
    y_mcl = np.array([halfheight]*len(x_mcl))


    geo_dict = {
        "periodics": {"ylower": yper_low, "yupper": yper_high},
        "span_z": span_z,
        "midpassagestreamLine": {"x_mpsl":x_mcl, "y_mpsl":y_mcl},
    }

    return geo_dict


def prep_geo(settings_yml):
    settings = yaml_dict_read(settings_yml)
    casepath = os.path.abspath(os.path.dirname(settings_yml))

    if "openfoam_channel" in settings["case_settings"]["case_type"]:
        geo_dict = prep_channelflow_geometry(settings_yml)

    geo_filename = "geometry.pkl"
    write_pickle(os.path.join(casepath, casedirs["data"], geo_filename), geo_dict)

    return geo_dict
