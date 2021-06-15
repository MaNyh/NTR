import matplotlib.pyplot as plt
import os
import csv

import numpy as np



def get_timeseries(file):
    with open(file, newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter='\t', quotechar='|')
        data = []
        for row in spamreader:
            data.append(row)
        data = data[5:]
        times = [float(i[0]) for i in data]
        massflow = [float(i[1]) for i in data]
    return np.array(times), np.array(massflow)


def plot():
    inlet_maindir = os.path.join("postProcessing", 'flowRatePatch(name=INLET)')
    inlet_massflow_dirs = os.listdir(inlet_maindir)
    inlet_dir = os.listdir(os.path.join(inlet_maindir, inlet_massflow_dirs[-1]))
    inlet_file = os.path.join(inlet_maindir, inlet_massflow_dirs[-1], inlet_dir[-1])

    outlet_maindir = os.path.join("postProcessing", 'flowRatePatch(name=OUTLET)')
    outlet_massflow_dirs = os.listdir(outlet_maindir)
    outlet_dir = os.listdir(os.path.join(outlet_maindir, outlet_massflow_dirs[-1]))
    outlet_file = os.path.join(outlet_maindir, outlet_massflow_dirs[-1], outlet_dir[-1])

    inlet_times, inlet_massflows = get_timeseries(inlet_file)
    outlet_times, outlet_massflows = get_timeseries(outlet_file)

    inlet_massflows = -inlet_massflows

    slicer = 1
    fig, axs = plt.subplots(3, 1)
    inlet_x = inlet_times[::slicer]
    inlet_y = inlet_massflows[::slicer]
    axs[0].plot(inlet_x, inlet_y, label="inlet")
    outlet_x = outlet_times[::slicer]
    outlet_y = outlet_massflows[::slicer]
    #axs[0].set_xlabel('time')
    axs[0].set_ylabel('massflow')
    axs[0].plot(outlet_x, outlet_y, label="outlet")
    axs[0].legend()

    residuals_maindir = os.path.join("postProcessing", "residuals")
    residuals_subdirs = os.listdir(residuals_maindir)
    res_dir = os.listdir(os.path.join(residuals_maindir, residuals_subdirs[-1]))
    res_file = os.path.join(residuals_maindir, residuals_subdirs[-1], res_dir[-1])

    res_times, residuals = get_timeseries(res_file)
    axs[1].plot(res_times, residuals,label="residuals")

    axs[1].set_xlabel('time')
    axs[1].set_ylabel('residuals')
    axs[1].legend()

    ypblade_maindir = os.path.join("postProcessing", "yPlusBlade")
    ypblade_subdirs = os.listdir(ypblade_maindir)
    ypb_dir = os.listdir(os.path.join(ypblade_maindir, ypblade_subdirs[-1]))
    ypb_file = os.path.join(ypblade_maindir, ypblade_subdirs[-1], ypb_dir[-1])

    ypb_times, yplusblade = get_timeseries(ypb_file)
    axs[2].plot(ypb_times, yplusblade,label="yplus blade")

    axs[2].set_xlabel('time')
    axs[2].set_ylabel('yPlus BLADE')
    axs[2].legend()




    plt.show()


plot()
