from NTR.postprocessing.openfoam.probe_monitoring import show_monitors
from NTR.postprocessing.openfoam.loginterpreter import logfilestats


show_monitors("case_settings.yml")
print(logfilestats("case_settings.yml"))
