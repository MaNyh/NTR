import os
import numpy as np
import matplotlib.pyplot as plt

from NTR.utils.filehandling import yaml_dict_read, read_csv
from NTR.database.case_dirstructure import casedirs
from NTR.utils.mathfunctions import vecAbs
from NTR.postprocessing.openfoam.loginterpreter import logfilestats


def show_monitors(casesettings_yml):
    settings = yaml_dict_read(casesettings_yml)
    if "INOUT_FIELDAVE_PROBING" in settings["simcase_settings"]["options"].keys():
        averagevaluesinlet(casesettings_yml)
        averagevaluesoutlet(casesettings_yml)
    if "INOUT_VELOCITY_PROBING" in settings["simcase_settings"]["options"]["INOUT_VELOCITY_PROBING"].keys():
        massflowoutlet(casesettings_yml)
        massflowinlet(casesettings_yml)
    if "XSCLICE_PROBING" in settings["simcase_settings"]["options"]["XSCLICE_PROBING"].keys():
        xslices(casesettings_yml)
    #if "logfile" in settings["monitoring"].keys():
    logfilestats("case_settings.yml")


def averagevaluesinlet(casesettings_yml):
    datname = "surfaceFieldValue"
    casepath = os.path.abspath(os.path.dirname(casesettings_yml))
    monitorpath = os.path.join(casepath,casedirs["solution"],"postProcessing","AverValuesInlet")
    make_averagevaluesplot(casepath, datname, monitorpath)


def averagevaluesoutlet(casesettings_yml):
    datname = "surfaceFieldValue"
    casepath = os.path.abspath(os.path.dirname(casesettings_yml))
    monitorpath = os.path.join(casepath,casedirs["solution"],"postProcessing","AverValuesOutlet")
    make_averagevaluesplot(casepath, datname, monitorpath)


def massflowoutlet(casesettings_yml):
    datname = "surfaceFieldValue"
    casepath = os.path.abspath(os.path.dirname(casesettings_yml))
    monitorpath = os.path.join(casepath,casedirs["solution"],"postProcessing","MassflowOutlet")
    make_averagevaluesplot(casepath, datname, monitorpath)


def massflowinlet(casesettings_yml):
    datname = "surfaceFieldValue"
    casepath = os.path.abspath(os.path.dirname(casesettings_yml))
    monitorpath = os.path.join(casepath,casedirs["solution"],"postProcessing","MassflowInlet")

    make_averagevaluesplot(casepath, datname, monitorpath)


def xslices(casesettings_yml):
    settings = yaml_dict_read(casesettings_yml)
    datnames = ["p", "U"]
    casepath = os.path.abspath(os.path.dirname(casesettings_yml))
    monitorpath = os.path.join(casepath,casedirs["solution"],"postProcessing","Probes_XSlices")

    nop = settings["simcase_optiondef"]["XSCLICE_PROBING"]["args"]["nop"]

    timeseries = {}

    for datname in datnames:
        dirlist = os.listdir(os.path.join(casepath, monitorpath))
        dirlistasfloats = [float(i) for i in dirlist]
        dirlistasfloats, dirlist = zip(*sorted(zip(dirlistasfloats, dirlist)))


        for timedirectory in dirlist:
            subdirlist = os.listdir(os.path.join(casepath, monitorpath, timedirectory))
            for item in subdirlist:
                if datname in item:
                    rawdata = read_csv(os.path.join(casepath, monitorpath, timedirectory, item))
                    if len(rawdata) <= nop + 2:
                        pass
                    else:
                        dats = rawdata[nop * 2 + 2:]

                        if datname not in timeseries.keys():
                            timeseries[datname] = {"time":[],
                                                   item:[]}

                        for row in dats:
                            if item=="U":
                                readablerow = row[0].split("                       ")
                                readablerow = list(filter(None, readablerow))
                                floatrow = readablerow[1:]
                                floatrow = [i.replace("(","").replace(")","").split(" ") for i in floatrow]
                                floatrow = [[float(y) for y in i] for i in floatrow]
                                floatrow = [vecAbs(np.array(i)) for i in floatrow]
                            else:
                                readablerow = row[0].split("       ")
                                readablerow = list(filter(None,readablerow))
                                floatrow = [float(i) for i in readablerow[1:]]

                            time = float(readablerow[0])
                            timeseries[datname][item].append(floatrow)
                            timeseries[datname]["time"].append(time)
    for item in subdirlist:
        timestamps = timeseries[item]["time"]
        bothslicevalues = np.asarray(timeseries[item][item])
        values = {"upstream": bothslicevalues[:,0:nop],
                  "downstream": bothslicevalues[:,nop:]}

        fig, axs = plt.subplots(1, 2)
        for idx, pos in enumerate(values.keys()):
            ys = values[pos]
            x = timestamps
            for idy in range(nop):
                y = ys[::,idy]
                axs[idx].plot(x, y, label=str(idy) + "-" + item)
            axs[idx].set_ylabel(pos + " " + item)
            for xvline in dirlistasfloats:
                axs[idx].axvline(x=xvline, linestyle=":", color="grey")
        fig.legend()
        fig.suptitle(monitorpath + " " + item)
        plt.show()
    return 0


def make_averagevaluesplot(casepath, datname, monitorpath):
    dirlist = os.listdir(os.path.join(casepath, monitorpath))
    dirlistasfloats = [float(i) for i in dirlist]
    dirlistasfloats, dirlist = zip(*sorted(zip(dirlistasfloats, dirlist)))
    timeseries = {"time": []}
    for timedirectory in dirlist:
        subdirlist = os.listdir(os.path.join(casepath, monitorpath, timedirectory))
        for item in subdirlist:
            if datname in item:

                rawdata = read_csv(os.path.join(casepath, monitorpath, timedirectory, item))
                if len(rawdata) <= 5:
                    pass

                else:
                    variables = rawdata[4][1:]

                    for var in variables:
                        if var not in timeseries.keys():
                            timeseries[var] = []

                    dat = np.asarray(rawdata[5:])

                    for idx, var in enumerate(timeseries.keys()):
                        if idx < dat.shape[1]:
                            newdat = dat[:, idx]
                        else:
                            newdat = ["0.0" for i in range(len(dat))]
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
    fig.suptitle(monitorpath)
    x = timeseries["time"]
    for idx, vars in enumerate(probe_variables):
        y = timeseries[vars]
        if len(probe_variables) > 1:
            x, y = zip(*sorted(zip(x, y)))
            axs[idx].plot(x, y, label=vars)
            axs[idx].set_ylabel(vars)
            for xvline in dirlistasfloats:
                axs[idx].axvline(x=xvline, linestyle=":", color="grey")
            axs[idx].legend()
        else:
            x, y = zip(*sorted(zip(x, y)))
            axs.plot(x, y, label=vars)
            axs.set_ylabel(vars)
            for xvline in dirlistasfloats:
                axs.axvline(x=xvline, linestyle=":", color="grey")
            axs.legend()

        # axs[idx].legend()
    plt.show()
