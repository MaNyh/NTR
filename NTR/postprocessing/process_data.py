import os

from NTR.utils.filehandling import yaml_dict_read
from NTR.postprocessing.process_from_settings import createProfileData_fromSettings

funcs = {
    "profile_data": createProfileData_fromSettings,
}


def postprocess(settings_yml):
    settings = yaml_dict_read(settings_yml)
    #casedir = os.path.dirname(settings_yml)
    casetype = settings["case_settings"]["type"]
    assert "post_settings" in settings.keys()
    funcs_to_call = settings["post_settings"]["call_funcs"]
    if casetype == "simulation":
        for f in funcs_to_call:
            funcs[f](settings_yml)
    elif casetype == "parastud":
        print(casetype)

