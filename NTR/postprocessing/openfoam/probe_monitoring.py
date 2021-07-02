import os
import numpy as np
import matplotlib.pyplot as plt


from NTR.utils.functions import yaml_dict_read, read_csv
from NTR.postprocessing.openfoam.loginterpreter import logfilestats

def show_monitors(casesettings_yml):
    settings = yaml_dict_read(casesettings_yml)
    if "inletoutletfieldave_probing" in settings["probing"]["probes"].keys():
        averagevaluesinlet(casesettings_yml)
        averagevaluesoutlet(casesettings_yml)
    if "inletoutletvelocity_probing" in settings["probing"]["probes"].keys():
        massflowoutlet(casesettings_yml)
        massflowinlet(casesettings_yml)
    if "logfile" in settings["monitoring"].keys():
        logfilestats("case_settings.yml")

def averagevaluesinlet(casesettings_yml):
    settings = yaml_dict_read(casesettings_yml)
    datname = "surfaceFieldValue"
    casepath = os.path.abspath(os.path.dirname(casesettings_yml))
    monitorpath = os.path.join(*settings["monitoring"]["averagevaluesinlet"].split("/"))

    make_averagevaluesplot(casepath, datname, monitorpath)

def averagevaluesoutlet(casesettings_yml):
    settings = yaml_dict_read(casesettings_yml)
    datname = "surfaceFieldValue"
    casepath = os.path.abspath(os.path.dirname(casesettings_yml))
    monitorpath = os.path.join(*settings["monitoring"]["averagevaluesoutlet"].split("/"))

    make_averagevaluesplot(casepath, datname, monitorpath)

def massflowoutlet(casesettings_yml):
    settings = yaml_dict_read(casesettings_yml)
    datname = "surfaceFieldValue"
    casepath = os.path.abspath(os.path.dirname(casesettings_yml))
    monitorpath = os.path.join(*settings["monitoring"]["massflowoutlet"].split("/"))

    make_averagevaluesplot(casepath, datname, monitorpath)

def massflowinlet(casesettings_yml):
    settings = yaml_dict_read(casesettings_yml)
    datname = "surfaceFieldValue"
    casepath = os.path.abspath(os.path.dirname(casesettings_yml))
    monitorpath = os.path.join(*settings["monitoring"]["massflowinlet"].split("/"))

    make_averagevaluesplot(casepath, datname, monitorpath)


def make_averagevaluesplot(casepath, datname, monitorpath):
    dirlist = os.listdir(os.path.join(casepath, monitorpath))
    timeseries = {"time": []}
    for timedirectory in dirlist:
        subdirlist = os.listdir(os.path.join(casepath, monitorpath, timedirectory))
        for item in subdirlist:
            if datname in item:

                rawdata = read_csv(os.path.join(casepath, monitorpath, timedirectory, item))
                variables = rawdata[4][1:]

                for var in variables:
                    if var not in timeseries.keys():
                        timeseries[var] = []

                dat = np.asarray(rawdata[5:])

                for idx, var in enumerate(timeseries.keys()):
                    newdat = dat[:, idx]
                    for row in newdat:
                        if "(" in row:
                            row = row.replace("(", "")
                            row = row.replace(")", "")
                            row = row.split(" ")
                            row = [float(i) for i in row]
                            timeseries[var].append(row)
                        else:
                            timeseries[var].append(float(row))
    probe_variables = list(timeseries.keys())[1:]

    fig, axs = plt.subplots(len(probe_variables), 1)
    x = timeseries["time"]
    for idx, vars in enumerate(probe_variables):
        y = timeseries[vars]
        if len(probe_variables)>1:
            axs[idx].plot(x, y, label=vars, marker="x",linestyle = 'None')
            axs[idx].set_ylabel(vars)
            axs[idx].legend()
        else:
            axs.plot(x, y, label=vars, marker="x",linestyle = 'None')
            axs.set_ylabel(vars)
            axs.legend()
        plt.title(monitorpath)
        #axs[idx].legend()
    plt.show()
