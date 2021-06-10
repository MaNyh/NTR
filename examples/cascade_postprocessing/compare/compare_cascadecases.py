from NTR.utils.casereader import case_from_dict
from NTR.postprocessing.createProfileData import createProfileData

rans_case = case_from_dict("../gwk_verdichterkaskade_rans/gwk_verdichtercascade_RANS.yml")
createProfileData(rans_case)

les_case = case_from_dict("../gwk_verdichterkaskade_les/gwk_verdichtercascade_LES.yml")
createProfileData(les_case)
