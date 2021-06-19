from NTR.preprocessing.openfoam.create_cascadecase_les import create_cascadecase_les
from NTR.preprocessing.openfoam.create_probes import create_probe_dicts
from NTR.utils.functions import yaml_dict_read

import os

def create_case(case_yml, global_yml):

    settings = yaml_dict_read(case_yml)
    globalparas = yaml_dict_read(global_yml)

    mainpath = os.path.abspath(os.path.dirname(case_yml))

    case_settings = settings["case"]
    probe_settings = settings["probes"]

    for k, v in globalparas.items():

        probe_settings[k] = v

    create_cascadecase_les(case_settings, mainpath, globalparas)
    create_probe_dicts(probe_settings)
