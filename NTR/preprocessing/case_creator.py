from NTR.preprocessing.openfoam.create_cascadecase_les import create_cascadecase_les
from NTR.preprocessing.openfoam.create_probes import create_probe_dicts
from NTR.utils.functions import yaml_dict_read

import os

def create_case(case_yml):
    settings = yaml_dict_read(case_yml)
    mainpath = os.path.abspath(os.path.dirname(case_yml))
    simcase_settings = settings["case"]
    probe_settings = settings["probes"]
    create_cascadecase_les(simcase_settings, mainpath)
    create_probe_dicts(probe_settings)
