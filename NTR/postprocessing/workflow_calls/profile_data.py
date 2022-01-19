import os
import matplotlib.pyplot as plt
import copy

from NTR.utils.filehandling import write_pickle, read_pickle
from NTR.utils.mesh_handling.pyvista_utils import load_mesh
from NTR.postprocessing.turbo.createProfileData import createProfileData


def profile_data_workflow(input, output):
    # =============================================================================
    # Daten Einlesen
    # =============================================================================
    # profile_data = {}

    refmesh = load_mesh(input)

    midspan_z = (refmesh.bounds[5] - refmesh.bounds[4]) / 2
    alpha = 0.005
    post_slice_1_x = refmesh.bounds[0] + 1e-6
    post_slice_2_x = refmesh.bounds[1] - 1e-6
    kappa = 1.4
    As = 1.458e-06
    Ts = 110.4
    R_L = 287
    outflow = refmesh.slice(normal="x", origin=(post_slice_2_x, 0, 0)).compute_cell_sizes()
    outflow = outflow.point_data_to_cell_data()
    p_k = outflow["p"] * outflow["Area"] / (sum(outflow["Area"]))
    cp = 1004.5
    l = 0.13
    output_path = os.path.dirname(input)

    profile_data = createProfileData(refmesh, midspan_z, alpha, post_slice_1_x, post_slice_2_x,
                                     output_path,
                                     kappa, R_L,
                                     p_k, As, l, cp, Ts)



    write_pickle(output, profile_data)


def plot_profiledata(inputs, output):
    fig, axs = plt.subplots(1, 3)

    lines = {}
    for respath in inputs:
        fname = os.path.basename(respath)[:-17]

        name, yangle_raw, prod_raw = fname.split("_")
        yangle = int(yangle_raw) / 100
        prod = int(prod_raw) / 10

        line_name = name + "_" + prod_raw

        if line_name not in lines.keys():
            lines[line_name] = {"mred": [], "pi": [], "eta_is": [], "lift_coefficient": []}
        res = read_pickle(respath)
        lines[line_name]["mred"].append(res["mred"])
        lines[line_name]["pi"].append(res["inte_p2"] / res["inte_p1"])
        lines[line_name]["eta_is"].append(res["inte_tot_p2"]-res["inte_tot_p1"])
        lines[line_name]["lift_coefficient"].append(res["lift_coefficient"])

    for name, line in lines.items():
        axs[0].plot(line["mred"], line["pi"], label=name)
        axs[1].plot(line["mred"], line["tot_p_loss"], label=name)
        axs[2].plot(line["beta1"]-line["beta1"], line["lift_coefficient"], label=name)
    axs[0].set_xlabel('mred')
    axs[0].set_ylabel('Pi')
    axs[1].set_xlabel('mred')
    axs[1].set_ylabel('tot_p_loss')
    axs[2].set_xlabel('mred')
    axs[2].set_ylabel('lift')
    # plt.show()
    plt.legend()
    plt.savefig(output)
    plt.close()


def plot_profilepressure(inputs, output):
    fig, axs = plt.subplots(1, 2)

    curves = {"angle":[],"ps_curve":[],"ss_curve":[]}

    curves_sets = {"reference": copy.deepcopy(curves), "ogrid": copy.deepcopy(curves), "wake": copy.deepcopy(curves)}
    curve_names = {"ogrid":0,"wake":1,"reference":-1}

    for respath in inputs:
        fname = os.path.basename(respath)[:-17]
        res = read_pickle(respath)

        name, yangle_raw, prod_raw = fname.split("_")
        if name not in curves.keys():
            curves[name] = {"ps":{},"ss":{}}
        casename = name + "_" + prod_raw
        prod = int(prod_raw)/100

        curves_sets[name]["angle"].append(yangle_raw)
        curves_sets[name]["ps_curve"].append([res["profileData"]["x_zu_l_ax_ps"], res["profileData"]["cp_ps"]])
        curves_sets[name]["ss_curve"].append([res["profileData"]["x_zu_l_ax_ss"], res["profileData"]["cp_ss"]])

    for curvetype, curvevals in curves_sets.items():
        ax_id = curve_names[curvetype]
        for idx in range(len(curvevals)):

            angle = curvevals["angle"][idx]
            ps_vals = curvevals["ps_curve"][idx]
            ss_vals = curvevals["ss_curve"][idx]

            if ax_id == -1:
                axs[0].plot(ps_vals[0],ps_vals[1], linestyle="dashed", color=(1, 1, 0), label=casename)
                axs[0].plot(ss_vals[0],ss_vals[1], linestyle="dashed", color=(1, 0, 0), label=casename)
                axs[1].plot(ps_vals[0],ps_vals[1], linestyle="dashed", color=(1, 1, 0), label=casename)
                axs[1].plot(ss_vals[0],ss_vals[1], linestyle="dashed", color=(1, 0, 0), label=casename)
            else:
                axs[ax_id].plot(ps_vals[0],ps_vals[1], linestyle="solid", color=(1,1,0), label=casename)
                axs[ax_id].plot(ss_vals[0],ss_vals[1], linestyle="solid", color=(1,0,0), label=casename)

   # plt.legend()
    plt.savefig(output)
    plt.close()

