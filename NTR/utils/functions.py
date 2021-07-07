import os
import yaml
import pickle
import csv


import NTR
from NTR.preprocessing.create_geom import create_geometry


def yaml_dict_read(yml_file):
    args_from_yaml = {}

    with open(yml_file, "r", newline='') as Fobj:
        document = yaml.load_all(Fobj, Loader=yaml.FullLoader)
        for settings in document:
            for key, value in settings.items():
                args_from_yaml[key] = value
    return args_from_yaml

def read_csv(csv_filepath):
    with open(csv_filepath,"r", newline='') as csvobj:
        spamreader = csv.reader(csvobj, delimiter='\t', quotechar='|')
        data = []
        for row in spamreader:
            data.append(row)
    return data

def write_igg_config(file, args):
    with open(file, "wb") as Fobj:
        pickle.dump(args, Fobj, protocol=0)


def run_igg_meshfuncs(settings_yaml):


    case_path = os.path.abspath(os.path.dirname(settings_yaml))
    settings = yaml_dict_read(settings_yaml)
    meshpath = os.path.join(case_path, "01_Meshing")
    if not os.path.isdir(meshpath):
        os.mkdir(meshpath)
    print(os.path.abspath(case_path))

    print("create_geometry")
    ptstxtfile = os.path.join(os.path.abspath(case_path), settings["geom"]["ptcloud_profile"])
    create_geometry(ptstxtfile,
                    settings["geometry"]["beta_meta_01"],
                    settings["geometry"]["beta_meta_02"],
                    settings["geom"]["x_inlet"],
                    settings["geom"]["x_outlet"],
                    settings["geometry"]["pitch"],
                    settings["geom"]["ptcloud_profile_unit"],
                    settings["geom"]["shift_domain"])

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

    write_igg_config(args_dict_path, args)
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


def readtxtfile(path_to_file):
    basepath = os.path.abspath(os.path.dirname(__file__))
    with open(os.path.join(basepath,path_to_file), "r") as fobj:
        content = fobj.readlines()
    return content
