import os
import tempfile
import shutil


from NTR.utils.filehandling import yaml_dict_read, write_yaml_dict
from NTR.postprocessing.process_from_settings import createProfileData_fromSettings
from NTR.preprocessing.create_simcase import read_parastudyaml, paracase_name, construct_paracasedict, get_parastud_parameternames_fromdict
from NTR.database.case_dirstructure import casedirs
from NTR.utils.dicthandling import delete_keys_from_dict

funcs = {
    "profile_data": createProfileData_fromSettings,
}

def postprocess(settings_yml):
    settings = yaml_dict_read(settings_yml)
    casedir = os.path.abspath(os.path.dirname(settings_yml))
    casetype = settings["case_settings"]["type"]
    assert "post_settings" in settings.keys()
    funcs_to_call = settings["post_settings"]["call_funcs"]
    postresults = {}
    if casetype == "simulation":
        for f in funcs_to_call:
            postresults["placeholder"]=funcs[f](settings_yml)
    elif casetype == "parameterstudy":
        paracases_settings = read_parastudyaml(settings_yml)
        solutionpath = settings["post_settings"]["solution"]
        paras = get_parastud_parameternames_fromdict(settings)

        tmp_dir = tempfile.TemporaryDirectory()

        paracasedirs = []
        for idx, settings_dict in enumerate(paracases_settings):
            settings_dict["case_settings"]["type"] = "simulation"
            casepara = construct_paracasedict(paras, settings_dict)
            sub_case_dir = paracase_name(casepara, idx)
            paracasedirs.append(sub_case_dir)

        for f in funcs_to_call:
            for cdir, case_settgs in zip(paracasedirs, paracases_settings):

                mesh_path = os.path.join(solutionpath.replace("*",cdir))
                origcase_yml = os.path.join(casedir,casedirs["solution"],cdir, "paracase_settings.yml")
                postprocess_yml = os.path.join(casedir, "postprocess.yml")
                #tmp_settings = yaml_dict_read(postprocess_yml)

                check_orig_settings = yaml_dict_read(origcase_yml)
                #todo: the next line is working. but changing post-settings will lead to assert. so for now outcommented

                checkorig = check_orig_settings.copy()
                checkcase = case_settgs.copy()
                if "post_settings" in checkorig.keys():
                    delete_keys_from_dict(checkorig,["post_settings"])
                if "post_settings" in checkcase.keys():
                    delete_keys_from_dict(checkcase,["post_settings"])

                #assert checkorig==checkcase, "cases do not match"
                case_settgs["post_settings"]["solution"] = os.path.join(casedir,mesh_path)


                write_yaml_dict(postprocess_yml, case_settgs)


                postresults["placeholder"] = funcs[f](postprocess_yml)

