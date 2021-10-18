import pyvista as pv
import imageio
import os

from NTR.database.case_dirstructure import casedirs
from NTR.utils.filehandling import yaml_dict_read


def create(path_to_yaml_dict):
    settings = yaml_dict_read(path_to_yaml_dict)
    casepath = os.path.abspath(os.path.dirname(path_to_yaml_dict))

    val_dict = {"yplus":,}

    volmesh
    wallmesh

#calc yplus

#calc AVE massflow
#calc AVE massflow diff