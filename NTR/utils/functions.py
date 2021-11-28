import os
import pickle

import NTR
from NTR.utils.filehandling import yaml_dict_read, write_pickle_protocolzero

from NTR.preprocessing.options.openfoam import openFoam_createProbesInletOutlet
from NTR.preprocessing.options.openfoam import openFoam_create_inletoutletave_probe_dict
from NTR.preprocessing.options.openfoam import openFoam_createXSliceProbes
from NTR.preprocessing.options.openfoam import openFoam_createProbesProfileDict
from NTR.preprocessing.options.openfoam import openFoam_createProbesStreamlineDict
from NTR.preprocessing.options.openfoam import openFoam_create_vk_stagflow_probes
from NTR.preprocessing.options.openfoam import openfoam_createSlicingDict

def run_igg_meshfuncs(settings_yaml):

    case_path = os.path.abspath(os.path.dirname(settings_yaml))
    settings = yaml_dict_read(settings_yaml)
    externals_settings = yaml_dict_read(os.path.join(os.path.dirname(__file__),"externals","externals_settings.yml"))
    meshpath = os.path.join(case_path, "01_Meshing")
    if not os.path.isdir(meshpath):
        os.mkdir(meshpath)
    print(os.path.abspath(case_path))

    print("create_mesh")
    cwd = os.getcwd()
    os.chdir(externals_settings["igg"]["install_directory"])
    igg_exe = externals_settings["igg"]["executable"]

    ntrpath = os.path.dirname(os.path.abspath(NTR.__file__))

    script_path = os.path.join(ntrpath, "utils", "externals", "numeca_igg", "igg_cascade_meshcreator.py")
    args_dict_path = os.path.join(ntrpath, "utils", "externals", "numeca_igg", externals_settings["igg"]["argument_pickle_dict"])

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


def all_equal(iterable):
    return iterable.count(iterable[0]) == len(iterable)

def func_by_name(val):
    funs = {"openFoam_createProbesInletOutlet":openFoam_createProbesInletOutlet,
            "openFoam_create_inletoutletave_probe_dict":openFoam_create_inletoutletave_probe_dict,
            "openFoam_createXSliceProbes":openFoam_createXSliceProbes,
            "openFoam_createProbesProfileDict":openFoam_createProbesProfileDict,
            "openFoam_createProbesStreamlineDict":openFoam_createProbesStreamlineDict,
            "openFoam_create_vk_stagflow_probes":openFoam_create_vk_stagflow_probes,
            "openfoam_createSlicingDict":openfoam_createSlicingDict
            }
    fun = funs[val]
    return fun
