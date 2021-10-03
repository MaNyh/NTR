import numpy as np


from NTR.utils.geom_functions.profileparas import extract_profile_paras

pointfile = "profile_pointcloud.txt"
profilepoints = np.loadtxt(pointfile)

ss_poly, ps_poly, centerline, beta_01, beta_02 = extract_profile_paras(pointfile)

