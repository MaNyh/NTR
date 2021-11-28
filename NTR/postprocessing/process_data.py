import os
import warnings
from tqdm import tqdm

from NTR.utils.filehandling import yaml_dict_read, write_yaml_dict, read_pickle, write_pickle
from NTR.postprocessing.process_from_settings import createProfileData_fromSettings
from NTR.preprocessing.create_simcase import read_parastudyaml, paracase_name, construct_paracasedict, \
    get_parastud_parameternames_fromdict
from NTR.database.case_dirstructure import casedirs
from NTR.utils.dicthandling import delete_keys_from_dict
from NTR.utils.mesh_handling.pyvista_utils import load_mesh, mesh_scalar_gradients

funcs = {
    "profile_data": createProfileData_fromSettings,
    "profile_loading":None,
}


def postprocess(settings_yml):
    settings = yaml_dict_read(settings_yml)
    casedir = os.path.abspath(os.path.dirname(settings_yml))
    casetype = settings["case_settings"]["type"]
    assert "post_settings" in settings.keys()
    funcs_to_call = settings["post_settings"]["call_funcs"]

    postdat_path = os.path.join(casedir,casedirs["data"],"postprocessed.pkl")

    if os.path.isfile(postdat_path):
        postresults = read_pickle(postdat_path)
    else:
        postresults = {}

    if casetype == "simulation":
        print("postprocessing " + casedir)

        for f in funcs_to_call:
            volmesh = settings["post_settings"]["solution"]
            postprocess_func(casedir, f, volmesh, postdat_path, settings_yml, postresults)

    elif casetype == "parameterstudy":
        paracases_settings = read_parastudyaml(settings_yml)
        solutionpath = settings["post_settings"]["solution"]
        paras = get_parastud_parameternames_fromdict(settings)

        paracasedirs = []
        for idx, settings_dict in enumerate(paracases_settings):
            settings_dict["case_settings"]["type"] = "simulation"
            casepara = construct_paracasedict(paras, settings_dict)
            sub_case_dir = paracase_name(casepara, idx)
            paracasedirs.append(sub_case_dir)

        with tqdm(total=len(paracasedirs)) as pbar:
            for cdir, case_settgs in zip(paracasedirs, paracases_settings):
                pbar.set_description(cdir)
                for f in funcs_to_call:

                    mesh_path = os.path.join(solutionpath.replace("*", cdir))
                    origcase_yml = os.path.join(casedir, casedirs["solution"], cdir, "paracase_settings.yml")
                    postprocess_yml = os.path.join(casedir, "postprocess.yml")

                    check_orig_settings = yaml_dict_read(origcase_yml)

                    checkorig = check_orig_settings.copy()
                    checkcase = case_settgs.copy()
                    if "post_settings" in checkorig.keys():
                        delete_keys_from_dict(checkorig, ["post_settings"])
                    if "post_settings" in checkcase.keys():
                        delete_keys_from_dict(checkcase, ["post_settings"])

                    if not (checkorig == checkcase):
                        warnings.warn("case-settings from central defined case and the solution does not match!")
                    case_settgs["post_settings"]["solution"] = os.path.join(casedir, mesh_path)

                    write_yaml_dict(postprocess_yml, case_settgs)

                    postprocess_func(cdir, f, mesh_path, postdat_path, postprocess_yml, postresults)

                pbar.update(1)


def postprocess_func(cdir, f, mesh_path, postdat_path, postprocess_yml, postresults):
    mesh_md_fsize = os.path.getsize(mesh_path)
    mesh_md_ftime = os.path.getmtime(mesh_path)
    if cdir in postresults.keys() and f in postresults[cdir].keys() and \
        postresults[cdir][f]["mesh_md"]["fsize"] == mesh_md_fsize and \
        postresults[cdir][f]["mesh_md"]["ftime"] == mesh_md_ftime:
        print("skipping...")
    else:
        postresults[cdir] = {}
        mesh = load_mesh(mesh_path)
        mesh = mesh_scalar_gradients(mesh, "U")
        try:
            postresults[cdir][f] = funcs[f](mesh, postprocess_yml)
            postresults[cdir][f]["mesh_md"] = {"fsize": mesh_md_fsize,
                                               "ftime": mesh_md_ftime,
                                               }
            write_pickle(postdat_path, postresults)
        except:
            postresults[cdir][f] = -1

        del mesh
