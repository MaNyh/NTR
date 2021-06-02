import os
import yaml
import pickle

import NTR
from NTR.utils.create_geom import create

def yaml_dict_read(yml_file):
    args_from_yaml = {}

    with open(yml_file, "r") as Fobj:
        document = yaml.load_all(Fobj, Loader=yaml.FullLoader)
        for settings in document:
            for key, value in settings.items():
                args_from_yaml[key] = value
    return args_from_yaml

def write_igg_config(file, args):
    with open(file, "wb") as Fobj:
        pickle.dump(args, Fobj, protocol=0)


def run_igg_meshfuncs(case_path):
    settings = yaml_dict_read(os.path.join(case_path, "settings.yml"))
    print(os.path.abspath(case_path))

    print("create_geometry")
    ptstxtfile = os.path.join(os.path.abspath(case_path),settings["geom"]["ptcloud_profile"])
    create(ptstxtfile,
           settings["geom"]["beta_meta_01"],
           settings["geom"]["beta_meta_02"],
           settings["geom"]["x_inlet"],
           settings["geom"]["x_outlet"],
           settings["geom"]["pitch"],
           settings["geom"]["ptcloud_profile_unit"] )


    print("create_mesh")
    cwd = os.getcwd()
    os.chdir(settings["igg"]["install_directory"])
    igg_exe = settings["igg"]["executable"]

    ntrpath = os.path.dirname(os.path.abspath(NTR.__file__))

    script_path = os.path.join(ntrpath, "utils", "externals", "numeca_igg", "igg_cascade_meshcreator.py")
    args_dict_path = os.path.join(ntrpath, "utils", "externals", "numeca_igg", settings["igg"]["argument_pickle_dict"])

    point_cloud_path = os.path.join(case_path, "geom.dat")

    args = {}

    args["pointcloudfile"] = point_cloud_path
    args["add_path"] = ntrpath
    args["case_path"] = case_path
    for i in settings["mesh"]:
        args[i] = settings["mesh"][i]

    args["save_project"] = os.path.join(case_path, 'mesh.igg')
    args["save_fluent"] = os.path.join(case_path, "fluent.msh")

    write_igg_config(args_dict_path, args)
    os.system(igg_exe + " -batch -print -script " + script_path)
    os.chdir(cwd)



def read_pickle_args(path):
    filepath = os.path.join(path, "args.pkl")
    with open(filepath,"rb") as Fobj:
        dict = pickle.load(Fobj)
    return dict


def absVec(vec):
    return (vec[0]**2+vec[1]**2+vec[2]**2)**0.5


def absvec_array(array):
    return [absVec(vec) for vec in array]
