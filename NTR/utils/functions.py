import os
import pickle

import NTR
from NTR.preprocessing.create_geom import create_geometry
from NTR.utils.filehandling import yaml_dict_read, write_pickle_protocolzero


def run_create_geometry(settings_yaml):
    case_path = os.path.abspath(os.path.dirname(settings_yaml))
    settings = yaml_dict_read(settings_yaml)
    meshpath = os.path.join(case_path, "01_Meshing")
    if not os.path.isdir(meshpath):
        os.mkdir(meshpath)
    print(os.path.abspath(case_path))

    print("create_geometry")
    ptstxtfile = os.path.join(os.path.abspath(case_path), settings["geom"]["ptcloud_profile"])
    outpath = os.path.join(os.path.dirname(os.path.abspath(settings_yaml)),"00_Ressources","01_Geometry")
    create_geometry(ptstxtfile,
                    settings["geom"]["x_inlet"],
                    settings["geom"]["x_outlet"],
                    settings["geometry"]["pitch"],
                    settings["geom"]["ptcloud_profile_unit"],
                    settings["geom"]["shift_domain"],
                    settings["geometry"]["alpha"],
                    settings["geometry"]["midline_tolerance"],
                    settings["mesh"]["extrudeLength"],
                    outpath)
    return 0

def run_igg_meshfuncs(settings_yaml):

    case_path = os.path.abspath(os.path.dirname(settings_yaml))
    settings = yaml_dict_read(settings_yaml)
    meshpath = os.path.join(case_path, "01_Meshing")
    if not os.path.isdir(meshpath):
        os.mkdir(meshpath)
    print(os.path.abspath(case_path))

    print("create_mesh")
    cwd = os.getcwd()
    os.chdir(settings["igg"]["install_directory"])
    igg_exe = settings["igg"]["executable"]

    ntrpath = os.path.dirname(os.path.abspath(NTR.__file__))

    script_path = os.path.join(ntrpath, "utils", "externals", "numeca_igg", "igg_cascade_meshcreator.py")
    args_dict_path = os.path.join(ntrpath, "utils", "externals", "numeca_igg", settings["igg"]["argument_pickle_dict"])

    point_cloud_path = os.path.join(case_path,"01_Meshing", "geom.dat")

    args = {"pointcloudfile": point_cloud_path,
            "add_path": ntrpath,
            "case_path": case_path,
            "save_project": os.path.join(case_path, 'mesh.igg'),
            "save_fluent": os.path.join(case_path, "fluent.msh")}

    for i in settings["mesh"]:
        args[i] = settings["mesh"][i]

    write_pickle_protocolzero(args_dict_path, args)
    os.system(igg_exe + " -batch -print -script " + script_path)
    os.chdir(cwd)


def read_pickle_args(path):
    filepath = os.path.join(path, "args.pkl")
    with open(filepath, "rb") as Fobj:
        pargs = pickle.load(Fobj)
    return pargs


def absVec(vec):
    return (vec[0]**2+vec[1]**2+vec[2]**2)**0.5


def absvec_array(array):
    return [absVec(vec) for vec in array]


