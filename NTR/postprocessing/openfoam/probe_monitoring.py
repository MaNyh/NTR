import os
import numpy as np
import matplotlib.pyplot as plt
import re

from NTR.utils.filehandling import yaml_dict_read, read_csv
from NTR.database.case_dirstructure import casedirs
from NTR.utils.mathfunctions import vecAbs
from NTR.postprocessing.openfoam.loginterpreter import logfilestats
from NTR.postprocessing.openfoam.probe_reading import readprobes

def show_monitors(casesettings_yml):
    settings = yaml_dict_read(casesettings_yml)
    if "INOUT_FIELDAVE_PROBING" in settings["simcase_settings"]["options"].keys():
        averagevaluesinlet(casesettings_yml)
        averagevaluesoutlet(casesettings_yml)
    if "INOUT_VELOCITY_PROBING" in settings["simcase_settings"]["options"].keys():
        massflowoutlet(casesettings_yml)
        massflowinlet(casesettings_yml)
    if "XSCLICE_PROBING" in settings["simcase_settings"]["options"].keys():
        xslices(casesettings_yml)
    #if "logfile" in settings["monitoring"].keys():
    logfilestats("case_settings.yml")


def averagevaluesinlet(casesettings_yml):
    datname = "surfaceFieldValue"
    casepath = os.path.abspath(os.path.dirname(casesettings_yml))
    monitorpath = os.path.join(casepath,casedirs["solution"],"postProcessing","AverValuesInlet")
    probe_variables = {"areaAverage(U)":[0,2], "areaAverage(p)":[3,3], "areaAverage(T)":[4,4]}
    make_averagevaluesplot(casepath, datname, monitorpath,probe_variables)


def averagevaluesoutlet(casesettings_yml):
    datname = "surfaceFieldValue"
    casepath = os.path.abspath(os.path.dirname(casesettings_yml))
    monitorpath = os.path.join(casepath,casedirs["solution"],"postProcessing","AverValuesOutlet")
    probe_variables = {"areaAverage(U)":[0,2], "areaAverage(p)":[3,3], "areaAverage(T)":[4,4]}
    make_averagevaluesplot(casepath, datname, monitorpath,probe_variables)


def massflowoutlet(casesettings_yml):
    datname = "surfaceFieldValue"
    casepath = os.path.abspath(os.path.dirname(casesettings_yml))
    monitorpath = os.path.join(casepath,casedirs["solution"],"postProcessing","MassflowOutlet")
    probe_variables = {"phi":[0,0]}
    make_averagevaluesplot(casepath, datname, monitorpath,probe_variables)


def massflowinlet(casesettings_yml):
    datname = "surfaceFieldValue"
    casepath = os.path.abspath(os.path.dirname(casesettings_yml))
    monitorpath = os.path.join(casepath,casedirs["solution"],"postProcessing","MassflowInlet")
    probe_variables = {"phi":[0,0]}
    make_averagevaluesplot(casepath, datname, monitorpath,probe_variables)


def xslices(casesettings_yml):
    settings = yaml_dict_read(casesettings_yml)
    datnames = ["pMean", "UPrime2Mean", "UMean"]
    casepath = os.path.abspath(os.path.dirname(casesettings_yml))
    monitorpath = os.path.join(casepath, casedirs["solution"], "postProcessing", "Probes_XSlices")

    if not os.path.isdir(monitorpath):
        print("no data found for func xslices")
        return 0

    nop = settings["simcase_optiondef"]["XSCLICE_PROBING"]["args"]["nop"]

    timeseries = {}

    for datname in datnames:
        dirlist = os.listdir(os.path.join(casepath, monitorpath))
        dirlistasfloats = [float(i) for i in dirlist]
        dirlistasfloats, dirlist = zip(*sorted(zip(dirlistasfloats, dirlist)))

        vectorpattern = "\(\-{0,1}\d{1,}.\d{1,}\s\-{0,1}\d{1,}.\d{1,}\s\-{0,1}\d{1,}.\d{1,}"

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

                            if re.search(vectorpattern, row[0]):
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


def make_averagevaluesplot(casepath, datname, monitorpath,probe_variables):
    dirlist = os.listdir(os.path.join(casepath, monitorpath))
    dirlistasfloats = [float(i) for i in dirlist]
    dirlistasfloats, dirlist = zip(*sorted(zip(dirlistasfloats, dirlist)))
    timeseries = {"time": []}
    for varname in probe_variables.keys():
        timeseries[varname]=[]

    for timedirectory in dirlist:
        subdirlist = os.listdir(os.path.join(casepath, monitorpath, timedirectory))
        for item in subdirlist:
            if datname in item:
                timesteps,vals_raw = readprobes(os.path.join(casepath, monitorpath),timedirectory,"","surfaceFieldValue.dat")
                vals_asarray = np.array([np.array(*i) for i in vals_raw])
                if len(vals_raw) <= 0:
                    pass
                timeseries["time"]=[*timeseries["time"],*timesteps]
                for var in probe_variables.keys():
                    lower = probe_variables[var][0]
                    upper = probe_variables[var][1]+1
                    vals_chunk = vals_asarray[:,lower:upper]
                    if var == "areaAverage(U)":
                        vals_chunk=[vecAbs(i) for i in vals_chunk]
                    timeseries[var] = [*timeseries[var], *vals_chunk]

    fig, axs = plt.subplots(len(probe_variables), 1)
    fig.suptitle(monitorpath)
    x = timeseries["time"]
    for idx, vars in enumerate(probe_variables.keys()):
        y = timeseries[vars]
        if len(y) > 1:
            x, y = zip(*sorted(zip(x, y)))
            #todo: error here. need fix
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
