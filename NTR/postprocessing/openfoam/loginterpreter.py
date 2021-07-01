import re
import numpy as np
import os

from NTR.utils.functions import readtxtfile, yaml_dict_read

class timestep:
    def __init__(self,raw_lines):
        pattern_cputime_line = "ExecutionTime"
        pattern_executiontime = "ExecutionTime\s{1,}=\s{1,}\d{,}.\d{1,}\s{1,}s"
        pattern_clocktime = "ClockTime\s{1,}=\s{1,}\d{,}\s{,}s"

        pattern_courant_line = "Courant Number"
        pattern_courant_mean = "Courant Number mean:\s{,}\d+\.\d+"
        pattern_courant_max = "max:\s\d{1,}\.\d{1,}"

        pattern_simulationtime_line = "Time\s=\s\d{1,}.\d{1,}(e-\d{1,}){,}"
        pattern_simulationtime = "\d{1,}.\d{1,}(e\d{1,}){,}"
        self.executiontime = 0
        self.clocktime = 0

        self.courant_mean = 0
        self.courant_max = 0

        self.simulationtime = 0

        for line in raw_lines:
            cputime_line_match = re.search(pattern_cputime_line, line)
            if cputime_line_match:
                executiontime_match = re.search(pattern_executiontime, line).span()
                self.executiontime = float(line[executiontime_match[0]+16:executiontime_match[1]-1])
                clocktime_match = re.search(pattern_clocktime, line).span()
                self.clocktime = float(line[clocktime_match[0]+12:clocktime_match[1]-1])

            courant_line_match = re.search(pattern_courant_line, line)
            if courant_line_match:
                courant_mean_match = re.search(pattern_courant_mean, line).span()
                self.courant_mean = float(line[courant_mean_match[0]+21:courant_mean_match[1]])
                courant_max_match = re.search(pattern_courant_max, line).span()
                self.courant_max = float(line[courant_max_match[0]+4:courant_max_match[1]])

            simulationtime_line_match = re.search(pattern_simulationtime_line, line)
            if simulationtime_line_match:
                simulationtime_match = re.search(pattern_simulationtime, line).span()
                self.simulationtime = float(line[simulationtime_match[0]:simulationtime_match[1]])


        self.raw_lines = raw_lines


def logfilestats(settings_dict):
    settings = yaml_dict_read(settings_dict)
    filepath = settings["monitoring"]["logfile"]
    casepath = os.path.abspath(os.path.dirname(settings_dict))
    logfile_raw = readtxtfile(os.path.join(casepath,filepath))
    timestepobj_list = []
    timestep_lines = []
    timestep_counter = 0
    for line in logfile_raw:

        timestep_lines.append(line)
        if re.search("Time\s=\s\d{1,}.\d{1,}(e-\d{1,}){,}", line):
            timestepobj_list.append(timestep(timestep_lines))
            timestep_counter += 1
            timestep_lines = []
    timestepobj_list.append(timestep(timestep_lines))
    timestep_counter += 1


    timesteps = len(timestepobj_list) - 3
    clocktime = timestepobj_list[-3].clocktime
    courant_max = max([i.courant_max for i in timestepobj_list])
    courant = np.mean([i.courant_max for i in timestepobj_list])
    timesteptime_mean = timesteps / clocktime

    return timesteps, clocktime, courant_max, courant, timesteptime_mean

