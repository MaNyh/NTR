from NTR.preprocessing.openfoam.create_cascadecase_les import create_cascadecase_les
from NTR.preprocessing.openfoam.create_probes import create_probe_dicts


def create_case(case_yml, probe_yml):
    create_cascadecase_les(case_yml)
    create_probe_dicts(probe_yml)
