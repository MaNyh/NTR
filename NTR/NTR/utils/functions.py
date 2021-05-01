import os
import yaml
import pickle

from NTR.utils.externals.tecplot_functions import openTecplotFile
from NTR.utils.create_geom import create

def yaml_dict_read(yml_file):

    args_from_yaml= {}

    with open(yml_file, "r") as Fobj:
        document = yaml.load_all(Fobj,Loader=yaml.FullLoader)
        for settings in document:
            for key, value in settings.items():
                args_from_yaml[key]=value

    return args_from_yaml

def write_igg_config(file, args):
    with open(file, "wb") as Fobj:
        pickle.dump(args, Fobj, protocol=0)


def run_igg_meshfuncs():
    global args
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
        script_path = os.path.join(cwd, ".." , "utils","externals", "igg_cascade_meshcreator.py")
        args_dict_path = os.path.join(cwd,".." , "utils","externals",settings["igg"]["argument_pickle_dict"])

        point_cloud_path = os.path.join(cwd, "ressources", "geom.dat")

        args = {}
        args["pointcloudfile"] = point_cloud_path
        args["openTecplotFile"] = openTecplotFile

        args["yPerLowHGridBlockPitchStart"] = settings["mesh"]["yPerLowHGridBlockPitchStart"]
        args["yPerHighHGridBlockPitchStart"] = settings["mesh"]["yPerHighHGridBlockPitchStart"]
        args["vk_BlockStartFromChord"] = settings["mesh"]["vk_BlockStartFromChord"]
        args["hk_BlockStartFromChord"] = settings["mesh"]["hk_BlockStartFromChord"]
        args["factor"] = settings["mesh"]["factor"]
        args["delta_i"] = settings["mesh"]["delta_i"]
        args["cellwidthcoeff"] = settings["mesh"]["cellwidthcoeff"]
        args["first_cell_width"] = settings["mesh"]["first_cell_width"]
        args["exp_ratio"] = settings["mesh"]["exp_ratio"]
        args["layers"] = settings["mesh"]["layers"]
        args["extrudeLength"] = settings["mesh"]["extrudeLength"]
        args["extrudeNodes"] = settings["mesh"]["extrudeNodes"]

        write_igg_config(args_dict_path, args)
        os.system(igg_exe + " -batch -print -script " + script_path)
        #os.remove(args_dict_path)
        # os.chdir(cwd)
    else:
        print("skipping meshing")
