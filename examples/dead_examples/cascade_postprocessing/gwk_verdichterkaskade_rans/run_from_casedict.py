from NTR.utils.casereader import case_from_dict
from NTR.postprocessing.createProfileData import createProfileData

c = case_from_dict("gwk_verdichtercascade_RANS.yml")
createProfileData(c)
