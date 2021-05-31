from NTR.utils.geom_functions import GetProfileValuesMidspan

def createProfileData(case):
    GetProfileValuesMidspan(case.mesh_dict["fluid"])
