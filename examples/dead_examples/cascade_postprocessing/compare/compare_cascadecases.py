from NTR.utils.casereader import case_from_dict
from NTR.postprocessing.createProfileData import createProfileData
import matplotlib.pyplot as plt
import os

thispath = os.getcwd()

fig = plt.figure()
rans_case = case_from_dict("../gwk_verdichterkaskade_rans/gwk_verdichtercascade_RANS.yml")
x_zu_l_ax_ss, cp_ss, x_zu_l_ax_ps, cp_ps = createProfileData(rans_case)
plt.plot(x_zu_l_ax_ss, cp_ss)
plt.plot(x_zu_l_ax_ps, cp_ps)

os.chdir(thispath)

les_case = case_from_dict("../gwk_verdichterkaskade_les_probecreation/gwk_verdichtercascade_LES.yml")
x_zu_l_ax_ss, cp_ss, x_zu_l_ax_ps, cp_ps = createProfileData(les_case)
plt.plot(x_zu_l_ax_ss, cp_ss)
plt.plot(x_zu_l_ax_ps, cp_ps)
plt.show()
