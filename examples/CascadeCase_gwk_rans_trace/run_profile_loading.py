import matplotlib.pyplot as plt

from NTR.postprocessing.profile_loading import calc_loading_volmesh
from NTR.database.datasets.measurement_data_2016111_gwk_compressor_gilge.interpret_raw import read_gilgegwk
from NTR.postprocessing.profile_loading import calc_inflow_cp

psVals, ssVals = calc_loading_volmesh("case_settings.yml")

# todo: this must be read by a generic function from the volmesh
fake_p1 = 100449.54753
fake_p1t = 102148.13498

ps_xc_numerical = psVals["xc"]
ps_xc_pressure_numerical = psVals["Pressure"]
ps_xc_cp_numerical = calc_inflow_cp(ps_xc_pressure_numerical, fake_p1t, fake_p1)

ss_xc_numerical = ssVals["xc"]
ss_xc_pressure_numerical = ssVals["Pressure"]
ss_xc_cp_numerical = calc_inflow_cp(ss_xc_pressure_numerical, fake_p1t, fake_p1)

ss_xc, ss_cp, ps_xc, ps_cp = read_gilgegwk(verbose=False)

plt.figure()
plt.plot(ss_xc,ss_cp , label="ss from experiment")
plt.plot( ps_xc,ps_cp, label="ps from experiment")
plt.plot(ss_xc_numerical, ss_xc_cp_numerical, label="ss from numeric")
plt.plot(ps_xc_numerical, ps_xc_cp_numerical, label="ss from numeric")
plt.show()
