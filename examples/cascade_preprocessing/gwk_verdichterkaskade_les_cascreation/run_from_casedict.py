from NTR.utils.casereader import case_from_dict
from NTR.preprocessing.openfoam.create_cascadecase_les import create_case

case = case_from_dict("gwk_verdichtercascade_LES.yml")
create_case(case)
