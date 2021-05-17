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
    #global args
    settings = yaml_dict_read("ressources/settings.yml")
    if settings["geom"]["create_bool"]:
        print("create_geometry")

        create(settings["geom"]["ptcloud_profile"],
               settings["geom"]["beta_meta_01"],
               settings["geom"]["beta_meta_02"],
               settings["geom"]["x_inlet"],
               settings["geom"]["x_outlet"],
               settings["geom"]["pitch"], )
    else:
        print("skipping geometry")
    if settings["mesh"]["create_bool"]:
        print("create_mesh")
        cwd = os.getcwd()
        os.chdir(settings["igg"]["install_directory"])
        igg_exe = settings["igg"]["executable"]

        ntrpath = os.path.dirname(os.path.abspath(NTR.__file__))

        script_path = os.path.join(ntrpath, "utils", "externals", "igg_cascade_meshcreator.py")
        args_dict_path = os.path.join(ntrpath, "tmp", settings["igg"]["argument_pickle_dict"])

        point_cloud_path = os.path.join(ntrpath,"..", "examples", "ressources", "geom.dat")

        args = {}

        args["pointcloudfile"] = point_cloud_path
        args["add_path"] = ntrpath
        args["case_path"] = case_path

        args["yPerLowHGridBlockPitchStart"] = settings["mesh"]["yPerLowHGridBlockPitchStart"]
        args["yPerHighHGridBlockPitchStart"] = settings["mesh"]["yPerHighHGridBlockPitchStart"]
        args["vk_BlockStartFromChord"] = settings["mesh"]["vk_BlockStartFromChord"]
        args["hk_BlockStartFromChord"] = settings["mesh"]["hk_BlockStartFromChord"]
        args["factor"] = settings["mesh"]["factor"]
        args["ogrid_factor"] = settings["mesh"]["ogrid_factor"]
        args["delta_i"] = settings["mesh"]["delta_i"]
        args["cellwidthcoeff"] = settings["mesh"]["cellwidthcoeff"]
        args["first_cell_width"] = settings["mesh"]["first_cell_width"]
        args["exp_ratio"] = settings["mesh"]["exp_ratio"]
        args["layers"] = settings["mesh"]["layers"]
        args["extrudeLength"] = settings["mesh"]["extrudeLength"]
        args["extrudeNodes"] = settings["mesh"]["extrudeNodes"]
        args["streamline_nodedensity_factor"] = settings["mesh"]["streamline_nodedensity_factor"]
        args["shift_vk_block_xaxiscoeff"] = settings["mesh"]["shift_vk_block_xaxiscoeff"]
        args["shift_hk_block_xaxiscoeff"] = settings["mesh"]["shift_hk_block_xaxiscoeff"]

        args["smoothing"] = settings["mesh"]["smoothing"]

        args["save_project"] = os.path.join(case_path, 'mesh.igg')
        args["save_fluent"] = os.path.join(case_path, "fluent.msh")

        write_igg_config(args_dict_path, args)
        os.system(igg_exe + " -print -script " + script_path)
        os.chdir(cwd)
    else:
        print("skipping meshing")


def read_pickle_args(path):
    filepath = os.path.join(path,"args.pkl")
    with open(filepath,"rb") as Fobj:
        dict = pickle.load(Fobj)
    return dict
