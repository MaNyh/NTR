from NTR.preprocessing.openfoam.create_foamcase import create_cascadecase_les, create_case
from NTR.preprocessing.openfoam.create_probes import create_probe_dicts
from NTR.utils.functions import yaml_dict_read

import os

def create_case(case_yml):

    settings = yaml_dict_read(case_yml)

    mainpath = os.path.abspath(os.path.dirname(case_yml))

    create_case(settings, mainpath)
    create_cascadecase_les(settings, mainpath)
    create_probe_dicts(settings)
