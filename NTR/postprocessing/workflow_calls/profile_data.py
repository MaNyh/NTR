import os
import matplotlib.pyplot as plt
import copy
import numpy as np
import pyvista as pv

from NTR.utils.filehandling import write_pickle, read_pickle
from NTR.utils.mesh_handling.pyvista_utils import load_mesh
from NTR.postprocessing.turbo.createProfileData import createProfileData




def profile_data_workflow(input, output):
    # =============================================================================
    # Daten Einlesen
    # =============================================================================
    # profile_data = {}
    refmesh = load_mesh(input)
    profile_data = {}

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
            lines[line_name] = {"mred": [], "pi": [], "tot_p_loss": [], "lift_coefficient": [], "beta1": []}
        res = read_pickle(respath)
        lines[line_name]["mred"].append(res["mred"])
        lines[line_name]["pi"].append(res["inte_p2"] / res["inte_p1"])
        lines[line_name]["tot_p_loss"].append(res["inte_p_tot2"] - res["inte_p_tot1"])
        lines[line_name]["lift_coefficient"].append(res["lift_coefficient"])
        lines[line_name]["beta1"].append(res["beta1"])

    for linename, line in lines.items():
        axs[0].plot(line["mred"], line["pi"], label=linename)
        axs[1].plot(line["mred"], line["tot_p_loss"], label=linename)
        axs[2].plot(line["beta1"], line["lift_coefficient"], label=linename)
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

def plot_profilepressure_comp(input, output):

        fname = os.path.basename(input)[:-17]
        name, yangle_raw, prod_raw = fname.split("_")
        print(input)
        print(fname)
        result_dict = read_pickle(input)
        refpath = os.path.join(os.path.dirname(input), "reference_" + yangle_raw + "_10_profile_data.pkl")
        reference_dict = read_pickle(refpath)

        casename = name + "_" + prod_raw
        prod = int(prod_raw) / 100


        plt.figure()

        if name == "reference":
            plt.plot(reference_dict["profileData"]["x_zu_l_ax_ps"], reference_dict["profileData"]["cp_ps"],label=name+"_ps" )
            plt.plot(reference_dict["profileData"]["x_zu_l_ax_ss"], reference_dict["profileData"]["cp_ss"],label=name+"_ss" )

        else:

            plt.plot(reference_dict["profileData"]["x_zu_l_ax_ps"], reference_dict["profileData"]["cp_ps"],
                     linestyle="solid", color=(1, 1, 0), label="reference")
            plt.plot(reference_dict["profileData"]["x_zu_l_ax_ss"], reference_dict["profileData"]["cp_ss"],
                     linestyle="solid", color=(1, 1, 0), label="reference")

            plt.plot(result_dict["profileData"]["x_zu_l_ax_ps"], result_dict["profileData"]["cp_ps"], linestyle="dashed",
                     color=(1, 1, 0), label=casename+"_ps")
            plt.plot(result_dict["profileData"]["x_zu_l_ax_ss"], result_dict["profileData"]["cp_ss"], linestyle="dashed",
                     color=(1, 1, 0), label=casename+"_ss")

        plt.legend()

        plt.savefig(output)
        plt.close()


def plot_profilepressure_all(inputs, output):
    fig, axs = plt.subplots(1, 2)

    curves = {"angle": [], "ps_curve": [], "ss_curve": []}

    curves_sets = {"reference": copy.deepcopy(curves), "ogrid": copy.deepcopy(curves), "wake": copy.deepcopy(curves)}
    curve_names = {"ogrid": 0, "wake": 1, "reference": -1}

    for respath in inputs:
        fname = os.path.basename(respath)[:-17]
        result_dict = read_pickle(respath)

        name, yangle_raw, prod_raw = fname.split("_")
        if name not in curves.keys():
            curves[name] = {"ps": {}, "ss": {}}
        casename = name + "_" + prod_raw
        prod = int(prod_raw) / 100

        curves_sets[name]["angle"].append(yangle_raw)
        curves_sets[name]["ps_curve"].append(
            [result_dict["profileData"]["x_zu_l_ax_ps"], result_dict["profileData"]["cp_ps"]])
        curves_sets[name]["ss_curve"].append(
            [result_dict["profileData"]["x_zu_l_ax_ss"], result_dict["profileData"]["cp_ss"]])

    for curvetype, curvevals in curves_sets.items():
        ax_id = curve_names[curvetype]
        for idx in range(len(curvevals)):

            angle = curvevals["angle"][idx]
            ps_vals = curvevals["ps_curve"][idx]
            ss_vals = curvevals["ss_curve"][idx]

            if ax_id == -1:
                axs[0].plot(ps_vals[0], ps_vals[1], linestyle="dashed", color=(1, 1, 0), label=casename+"_ps")
                axs[0].plot(ss_vals[0], ss_vals[1], linestyle="dashed", color=(1, 0, 0), label=casename+"_ss")
                axs[1].plot(ps_vals[0], ps_vals[1], linestyle="dashed", color=(1, 1, 0), label=casename+"_ps")
                axs[1].plot(ss_vals[0], ss_vals[1], linestyle="dashed", color=(1, 0, 0), label=casename+"_ss")
            else:
                axs[ax_id].plot(ps_vals[0], ps_vals[1], linestyle="solid", color=(1, 1, 0), label=casename+"_ps")
                axs[ax_id].plot(ss_vals[0], ss_vals[1], linestyle="solid", color=(1, 0, 0), label=casename+"_ss")

    # plt.legend()
    plt.savefig(output)
    plt.close()


def plot_entropy_comp(input, output):
    casename = input.split("/")[1].split("_")[0]
    refcase = input.split("/")[1].replace(casename,"reference")[:-2]+"10"
    refpath = os.path.join(*input.split("/")[:-4],refcase,"output","cgns","TRACE.cgns")

    resultmesh = load_mesh(input)
    resultmesh.rotate_z(90)
    bounds_z = resultmesh.bounds[5] - resultmesh.bounds[4]
    midspanplane_result = resultmesh.slice(normal="z", origin=(0, 0, bounds_z / 2))


    pv.set_plot_theme("document")

    if casename=="reference":
        compute_entropy(midspanplane_result)
        p = pv.Plotter(off_screen=True)
        p.add_mesh(midspanplane_result,scalars="s")
        p.show(screenshot=output,cpos=(0,0,1),window_size=[4800,4800])
    else:
        shift_y = resultmesh.bounds[3] - resultmesh.bounds[2]
        refmesh = load_mesh(refpath)
        refmesh.rotate_z(90)
        midspanplane_reference = refmesh.slice(normal="z", origin=(0, 0, bounds_z / 2))
        midspanplane_reference.translate((0,shift_y,0))
        compute_entropy(midspanplane_result)
        compute_entropy(midspanplane_reference)
        p = pv.Plotter(off_screen=True)
        p.add_mesh(midspanplane_reference,scalars="s")
        p.add_mesh(midspanplane_result,scalars="s")
        p.show(screenshot=output,cpos=(0,0,1),window_size=[4800,4800])

def compute_entropy(mesh):
    cp = 1.0035

    R = 8.31446261815324
    s0 = 0
    T0 = 293.15
    p0 = 1.0135e5
    s = s0 + cp * np.log(mesh["T"]/T0)-R*np.log(mesh["p"]/p0)
    mesh["s"]=s


def plot_entropy_comp_diff(input, output):
    casename = input.split("/")[1].split("_")[0]
    refcase = input.split("/")[1].replace(casename,"reference")[:-2]+"10"
    refpath = os.path.join(*input.split("/")[:-4],refcase,"output","cgns","TRACE.cgns")

    resultmesh = load_mesh(input)
    resultmesh.rotate_z(90)
    bounds_z = resultmesh.bounds[5] - resultmesh.bounds[4]
    midspanplane_result = resultmesh.slice(normal="z", origin=(0, 0, bounds_z / 2))

    pv.set_plot_theme("document")

    if casename=="reference":
        compute_entropy(midspanplane_result)
        p = pv.Plotter(off_screen=True)
        p.add_mesh(midspanplane_result,scalars="s")
        p.show(screenshot=output,cpos=(0,0,1),window_size=[4800,4800])
    else:
        refmesh = load_mesh(refpath)
        refmesh.rotate_z(90)
        midspanplane_reference = refmesh.slice(normal="z", origin=(0, 0, bounds_z / 2))
        compute_entropy(midspanplane_result)
        compute_entropy(midspanplane_reference)
        midspanplane_result["sdiff"] =midspanplane_result["s"]-midspanplane_reference["s"]
        p = pv.Plotter(off_screen=True)
        p.add_mesh(midspanplane_result,scalars="sdiff")
        p.show(screenshot=output,cpos=(0,0,1),window_size=[4800,4800])


def plot_countours(input,output_U,output_p,output_T,output_rho):
    casename = input.split("/")[1].split("_")[0]
    refcase = input.split("/")[1].replace(casename,"reference")[:-2]+"10"
    refpath = os.path.join(*input.split("/")[:-4],refcase,"output","cgns","TRACE.cgns")

    resultmesh = load_mesh(input)
    resultmesh.rotate_z(90)
    bounds_z = resultmesh.bounds[5] - resultmesh.bounds[4]
    midspanplane_result = resultmesh.slice(normal="z", origin=(0, 0, bounds_z / 2))

    pv.set_plot_theme("document")

    if casename=="reference":
        p = pv.Plotter(off_screen=True)
        p.add_mesh(midspanplane_result,scalars="U")
        p.show(screenshot=output_U,cpos=(0,0,1),window_size=[4800,4800])
    else:
        shift_y = resultmesh.bounds[3] - resultmesh.bounds[2]
        refmesh = load_mesh(refpath)
        refmesh.rotate_z(90)
        midspanplane_reference = refmesh.slice(normal="z", origin=(0, 0, bounds_z / 2))
        midspanplane_reference.translate((0,shift_y,0))

        p = pv.Plotter(off_screen=True)
        p.add_mesh(midspanplane_reference,scalars="U")
        p.add_mesh(midspanplane_result,scalars="U")
        p.show(screenshot=output_U,cpos=(0,0,1),window_size=[4800,4800])


    if casename=="reference":
        p = pv.Plotter(off_screen=True)
        p.add_mesh(midspanplane_result,scalars="p")
        p.show(screenshot=output_p,cpos=(0,0,1),window_size=[4800,4800])
    else:
        shift_y = resultmesh.bounds[3] - resultmesh.bounds[2]
        refmesh = load_mesh(refpath)
        refmesh.rotate_z(90)
        midspanplane_reference = refmesh.slice(normal="z", origin=(0, 0, bounds_z / 2))
        midspanplane_reference.translate((0,shift_y,0))

        p = pv.Plotter(off_screen=True)
        p.add_mesh(midspanplane_reference,scalars="p")
        p.add_mesh(midspanplane_result,scalars="p")
        p.show(screenshot=output_p,cpos=(0,0,1),window_size=[4800,4800])

    if casename=="reference":
        p = pv.Plotter(off_screen=True)
        p.add_mesh(midspanplane_result,scalars="T")
        p.show(screenshot=output_T,cpos=(0,0,1),window_size=[4800,4800])
    else:
        shift_y = resultmesh.bounds[3] - resultmesh.bounds[2]
        refmesh = load_mesh(refpath)
        refmesh.rotate_z(90)
        midspanplane_reference = refmesh.slice(normal="z", origin=(0, 0, bounds_z / 2))
        midspanplane_reference.translate((0,shift_y,0))

        p = pv.Plotter(off_screen=True)
        p.add_mesh(midspanplane_reference,scalars="T")
        p.add_mesh(midspanplane_result,scalars="T")
        p.show(screenshot=output_T,cpos=(0,0,1),window_size=[4800,4800])


    if casename=="reference":
        p = pv.Plotter(off_screen=True)
        p.add_mesh(midspanplane_result,scalars="rho")
        p.show(screenshot=output_rho,cpos=(0,0,1),window_size=[4800,4800])
    else:
        shift_y = resultmesh.bounds[3] - resultmesh.bounds[2]
        refmesh = load_mesh(refpath)
        refmesh.rotate_z(90)
        midspanplane_reference = refmesh.slice(normal="z", origin=(0, 0, bounds_z / 2))
        midspanplane_reference.translate((0,shift_y,0))

        p = pv.Plotter(off_screen=True)
        p.add_mesh(midspanplane_reference,scalars="rho")
        p.add_mesh(midspanplane_result,scalars="rho")
        p.show(screenshot=output_rho,cpos=(0,0,1),window_size=[4800,4800])
