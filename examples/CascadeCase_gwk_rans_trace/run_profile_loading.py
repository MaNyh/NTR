from NTR.postprocessing.profile_loading import calc_loading_volmesh
from NTR.database.datasets.measurement_data_2016111_gwk_compressor_gilge.interpret_raw import read_gilgegwk

#psVals, ssVals = calc_loading_volmesh("case_settings.yml")
ss_x, ss_cp, ps_x, ps_cp = read_gilgegwk()
