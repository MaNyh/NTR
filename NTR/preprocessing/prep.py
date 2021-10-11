import pyvista as pv
import os

from NTR.utils.filehandling import yaml_dict_read, write_pickle
from NTR.database.case_dirstructure import casedirs

def prep_channelflow_geometry(settings_yml):
    settings=yaml_dict_read(settings_yml)

    length = float(settings["geometry"]["channel_length"])
    width = float(settings["geometry"]["yper_width"])
    halfheight = float(settings["geometry"]["yper_width"])

    yper_low = pv.Line((0,0,0),(length,0,0))
    yper_high = pv.Line((0,width,0),(length,width,0))
    span_z = halfheight

    geo_dict = {
                "periodics": {"ylower": yper_low, "yupper": yper_high},
                "span_z" : span_z,
                }

    return geo_dict


def prep_geo(settings_yml):
    settings = yaml_dict_read(settings_yml)
    casepath = os.path.abspath(os.path.dirname(settings_yml))

    if settings["case_settings"]["case_type"] == "openfoam_channel_les":
        geo_dict = prep_channelflow_geometry(settings_yml)

    geo_filename = "geometry.pkl"
    write_pickle(os.path.join(casepath, casedirs["data"], geo_filename), geo_dict)

    return geo_dict