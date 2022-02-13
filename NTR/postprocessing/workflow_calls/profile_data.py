import os
import matplotlib.pyplot as plt
import copy
import numpy as np
import pandas as pd
import pyvista as pv
import numpy as np
import matplotlib.pylab as pl
from matplotlib.colors import ListedColormap
from scipy.interpolate import griddata

from NTR.utils.filehandling import write_pickle, read_pickle
from NTR.utils.mesh_handling.pyvista_utils import load_mesh
from NTR.postprocessing.turbo.createProfileData import createProfileData
from NTR.utils.mesh_handling.pyvista_utils import mesh_scalar_gradients
from NTR.utils.mathfunctions import radToGrad,vecAngle, vecAbs, vecProjection


def areaAvePlane(mesh, val):
    array = mesh[val]
    areas = mesh["Area"]
    area_ave = sum((array.T * areas).T) / sum(areas)
    return area_ave


def massflowPlane(mesh):
    if not "Normals" in mesh.array_names:
        mesh = mesh.compute_normals()
    if not "Area" in mesh.array_names:
        mesh = mesh.compute_cell_sizes()
    mesh = mesh.point_data_to_cell_data()
    normals = mesh.cell_normals
    rhos = mesh["rho"]
    areas = mesh["Area"]
    velocities = mesh["U"]

    massflow = np.array(
        [vecAbs(vecProjection(velocities[i], normals[i])) for i in range(mesh.number_of_cells)]) ** 2 * rhos * areas

    return massflow


def massflowAvePlane(mesh, val):
    massflow = massflowPlane(mesh)

    mass_ave = sum(mesh[val] * massflow) / sum(massflow)
    return mass_ave


def rigvals(inlet, outlet, blade, output, p_ref, T_ref):
    path = os.path.normpath(os.path.dirname(inlet))
    # casename, yangle, prod = os.path.basename(inlet).split("-")
    inlet = load_mesh(inlet).extract_surface()
    outlet = load_mesh(outlet).extract_surface()
    blade = load_mesh(blade).extract_surface()

    if not "Normals" in inlet.array_names:
        inlet = inlet.compute_normals()
    if not "Area" in inlet.array_names:
        inlet = inlet.compute_cell_sizes()

    if not "Normals" in blade.array_names:
        blade = blade.compute_normals()
    if not "Area" in blade.array_names:
        blade = blade.compute_cell_sizes()
    blade = blade.point_data_to_cell_data()

    inlet = inlet.point_data_to_cell_data()
    inlet_massflow = sum(massflowPlane(inlet))

    if not "Normals" in outlet.array_names:
        outlet = outlet.compute_normals()
    if not "Area" in outlet.array_names:
        outlet = outlet.compute_cell_sizes()
    outlet = outlet.point_data_to_cell_data()

    outlet_massflow = sum(massflowPlane(outlet))

    p_out = areaAvePlane(outlet, "p")
    p_out_dyn = np.array(outlet["U"] ** 2).sum(axis=1) * outlet["rho"] / 2
    outlet["p_out_dyn"] = p_out_dyn
    p_out_dyn_ave = areaAvePlane(outlet, "p_out_dyn")

    p_in = areaAvePlane(inlet, "p")
    p_in_dyn = np.array(inlet["U"] ** 2).sum(axis=1) * inlet["rho"] / 2
    inlet["p_in_dyn"] = p_in_dyn
    p_in_dyn_ave = areaAvePlane(inlet, "p_in_dyn")

    inlet["u-normal"] = np.array([vecAbs(vecProjection(u, n)) for u, n in zip(inlet["U"], inlet["Normals"])])
    inlet["i-normal"] = inlet["u-normal"] * inlet["rho"]
    m_s = areaAvePlane(inlet, "i-normal")  # areaAvePlane(inflow,"U")*areaAvePlane(inflow,"rho")
    m_red_s = ((m_s * (p_ref / areaAvePlane(inlet, "p"))) * (areaAvePlane(inlet, "T") / T_ref)) ** .5
    p_tot_loss = -((p_out + p_out_dyn_ave) - (p_in + p_in_dyn_ave))

    force = blade["Area"] * blade["p"]
    blade["force_y"] = blade["Normals"][:, 1] * force
    blade["force_x"] = blade["Normals"][:, 0] * force

    inte_rho_in = massflowAvePlane(inlet, "rho")
    inte_U_in = vecAbs(areaAvePlane(inlet, "U"))

    liftcoefficient = abs(sum(blade["force_y"])) / (0.5 * inte_rho_in * inte_U_in ** 2 * sum(blade["Area"]))

    # blade["force_z"] = [vecProjection(np.array([0, 0, 1]), i) for i in blade["Normals"] * blade["force"][:, np.newaxis]]

    rigvals = {"p_out": p_out,
               "p_dyn_out": p_out_dyn,
               "m_in": inlet_massflow,
               "p_in": p_in,
               "p_dyn_in": p_in_dyn,
               "m_out": outlet_massflow,
               "m_red": m_red_s,
               "p_tot_loss": p_tot_loss,
               "blade_forcey": sum(blade["force_y"]),
               "blade_forcex": sum(blade["force_x"]),
               "lift_coefficient": liftcoefficient,
               }

    write_pickle(output, rigvals)
    # mesh.save(output)


def zslice_domain(mesh, output, rheight):
    """

    :param input: os.path to input-file
    :param output: os.path to output-file
    :param relative_heights: list of floats

    """

    refmesh = mesh_scalar_gradients(load_mesh(mesh), "U")
    bounds = refmesh.bounds
    zspan = np.sign((bounds[5] - bounds[4])) * (bounds[5] - bounds[4]) * float(rheight) + bounds[4]
    slice = refmesh.slice(normal="z", origin=(0, 0, zspan))
    slice = slice.compute_normals()
    slice.save(output)


def profile_data_workflow(input, output, rigvals, alpha, kappa, As, Ts, R_L, cp):
    # =============================================================================
    # Daten Einlesen
    # =============================================================================
    slice = load_mesh(input)
    rigval = read_pickle(rigvals)
    p_k = rigval["p_out"]

    post_slice_1_x = slice.bounds[0] + 1e-6
    post_slice_2_x = slice.bounds[1] - 1e-6

    output_path = os.path.dirname(input)

    profile_data = createProfileData(slice, alpha, post_slice_1_x, post_slice_2_x,
                                     output_path,
                                     kappa, R_L,
                                     p_k, As, cp, Ts)

    write_pickle(output, profile_data)


def plot_profiledata(inputs, output, beta_1):
    fig, axs = plt.subplots(1, 3, figsize=(15, 8))

    lines = {}
    for respath in inputs:
        fname = os.path.basename(respath)

        name, yangle_raw, prod_raw = fname.split("-")[0], fname.split("-")[1], fname.split("-")[2]
        yangle = int(yangle_raw) / 100
        ## prod = int(prod_raw) / 10

        line_name = name + "_" + prod_raw

        if line_name not in lines.keys():
            lines[line_name] = {"mred": [], "pi": [], "tot_p_loss": [], "lift_coefficient": [], "alpha": []}
        res = read_pickle(respath)
        lines[line_name]["mred"].append(res["m_red"])
        lines[line_name]["pi"].append(res["p_out"] / res["p_in"])
        lines[line_name]["tot_p_loss"].append(res["p_tot_loss"])
        lines[line_name]["lift_coefficient"].append(abs(res["lift_coefficient"]))
        lines[line_name]["alpha"].append(yangle - beta_1)

    for linename, line in lines.items():
        axs[0].scatter(line["mred"], line["pi"], label=linename)
        axs[1].scatter(line["alpha"], line["tot_p_loss"], label=linename)
        axs[2].scatter(line["alpha"], line["lift_coefficient"], label=linename)

    axs[0].set_xlabel(r'$\dot{m}_{red}$', fontsize=20.0)
    axs[0].set_ylabel(r'$\Pi_{stat}$', fontsize=20.0)
    axs[1].set_xlabel(r'$\alpha$', fontsize=20.0)
    axs[1].set_ylabel(r'$\Delta p_{tot,v}$', fontsize=20.0)
    axs[2].set_xlabel(r'$\alpha$', fontsize=20.0)
    axs[2].set_ylabel(r'$c_l$', fontsize=20.0)

    axs[0].grid()
    axs[1].grid()
    axs[2].grid()

    plt.legend()
    plt.tight_layout()
    plt.savefig(output)
    plt.close()


def plot_profilepressure_comp(input, output, highlim, lowlim):
    fname = os.path.basename(input)[:-17]
    name, yangle_raw, prod_raw, zslice = fname.split("-")
    result_dict = read_pickle(input)

    refpath = os.path.join(input.replace(name, "reference").replace(prod_raw, "10"))

    reference_dict = read_pickle(refpath)

    casename = name + "_" + prod_raw

    plt.figure()

    if name == "reference":
        plt.plot(reference_dict["profileData"]["x_zu_l_ax_ps"], reference_dict["profileData"]["cp_ps"],
                 label=name + "_ps", color="blue")
        plt.plot(reference_dict["profileData"]["x_zu_l_ax_ss"], reference_dict["profileData"]["cp_ss"],
                 label=name + "_ss", color="red")

    else:

        plt.plot(reference_dict["profileData"]["x_zu_l_ax_ps"], reference_dict["profileData"]["cp_ps"],
                 linestyle="solid", label="reference_ps", color="blue")
        plt.plot(reference_dict["profileData"]["x_zu_l_ax_ss"], reference_dict["profileData"]["cp_ss"],
                 linestyle="solid", label="reference_ss", color="red")

        plt.plot(result_dict["profileData"]["x_zu_l_ax_ps"], result_dict["profileData"]["cp_ps"], linestyle="dashed",
                 label=casename + "_ps", color="blue")
        plt.plot(result_dict["profileData"]["x_zu_l_ax_ss"], result_dict["profileData"]["cp_ss"], linestyle="dashed",
                 label=casename + "_ss", color="red")
    plt.grid(True, linestyle='-', linewidth=1)
    plt.legend()
    plt.ylim((lowlim, highlim))

    plt.title(input)
    plt.savefig(output)
    plt.close()


def plot_entropy_comp(input, output):
    casename = input.split("/")[1].split("_")[0]
    refcase = input.split("/")[1].replace(casename, "reference")[:-2] + "10"
    refpath = os.path.join(*input.split("/")[:-4], refcase, "output", "cgns", "TRACE.cgns")

    resultmesh = load_mesh(input)
    resultmesh.rotate_z(90)
    bounds_z = resultmesh.bounds[5] - resultmesh.bounds[4]

    pv.set_plot_theme("document")

    if casename == "reference":
        compute_entropy(resultmesh)
        p = pv.Plotter(off_screen=True)
        p.add_mesh(resultmesh, scalars="s")
        p.show(screenshot=output, cpos=(0, 0, 1), window_size=[4800, 4800], title=input)
    else:
        shift_y = resultmesh.bounds[3] - resultmesh.bounds[2]
        refmesh = load_mesh(refpath)
        refmesh.rotate_z(90)
        midspanplane_reference = refmesh.slice(normal="z", origin=(0, 0, bounds_z / 2))
        midspanplane_reference.translate((0, shift_y, 0))
        compute_entropy(midspanplane_result)
        compute_entropy(midspanplane_reference)
        p = pv.Plotter(off_screen=True)
        p.add_mesh(midspanplane_reference, scalars="s")
        p.add_mesh(midspanplane_result, scalars="s")
        p.show(screenshot=output, cpos=(0, 0, 1), window_size=[4800, 4800], title=input)


def compute_entropy(mesh, cp, R, Tref, pref):
    s0 = 0
    s = s0 + cp * np.log(mesh["T"] / Tref) - R * np.log(mesh["p"] / pref)
    mesh["s"] = s


def plot_entropy_comp_diff(input, output, cp, R, Tref, pref, s_max, s_min, sdiff_max, sdiff_min):
    input_split = os.path.basename(input).split("-")
    casename, yangle, prod, zslice = input_split[0], input_split[1], input_split[2], input_split[3]
    casename = input.split("/")[-1].split("-")[0]
    refpath = input.replace(casename, "reference").replace(prod, "10")

    resultmesh = load_mesh(input)
    resultmesh.rotate_z(90)

    pv.set_plot_theme("document")

    res = 4800
    title_size = int(0.02 * res)
    sargs = dict(
        title_font_size=title_size,
        label_font_size=int(0.016 * res),
        shadow=True,
        n_labels=3,
        italic=True,
        # fmt="%.1f",
        font_family="arial",
    )

    p = pv.Plotter(off_screen=True)
    p.add_title(input, font_size=title_size)
    if casename == "reference":
        compute_entropy(resultmesh, cp, R, Tref, pref)
        p.add_mesh(resultmesh, scalars="s", scalar_bar_args=sargs, cmap="coolwarm", clim=[s_min, s_max])
        p.show(screenshot=output, cpos=(0, 0, 1), window_size=[res, res])
    else:
        refmesh = load_mesh(refpath)
        refmesh.rotate_z(90)
        compute_entropy(resultmesh, cp, R, Tref, pref)
        compute_entropy(refmesh, cp, R, Tref, pref)
        resultmesh["sdiff"] = resultmesh["s"] - refmesh["s"]

        p.add_mesh(resultmesh, scalars="sdiff", scalar_bar_args=sargs, cmap="coolwarm", clim=[sdiff_min, sdiff_max])
        p.show(screenshot=output, cpos=(0, 0, 1), window_size=[res, res])


def plot_countours(input, output_U, output_p, output_T, output_rho,output_tke,
        vel_low,vel_high,tke_low,tke_high,
        temperature_low,temperature_high,rho_low,
        rho_high):
    directory = input.split("/")[0]
    casename = input.split("/")[1].split("-")[0]
    prod = input.split("/")[1].split("-")[2]
    refcase = input.split("/")[1].replace(casename, "reference").replace(prod, "10")
    refpath = os.path.join(directory, refcase)

    resultmesh = load_mesh(input)
    resultmesh.rotate_z(90)
    pv.set_plot_theme("document")

    res = 4800
    title_size = int(0.02 * res)
    sargs = dict(
        title_font_size=title_size,
        label_font_size=int(0.016 * res),
        shadow=True,
        n_labels=3,
        italic=True,
        # fmt="%.1f",
        font_family="arial",
    )
    if casename == "reference":
        p = pv.Plotter(off_screen=True)
        p.add_mesh(resultmesh, scalars="U", scalar_bar_args=sargs, cmap="coolwarm",clim=(vel_low,vel_high) )
        p.show(screenshot=output_U, cpos=(0, 0, 1), window_size=[4800, 4800])


    else:
        shift_y = resultmesh.bounds[3] - resultmesh.bounds[2]
        refmesh = load_mesh(refpath)

        refmesh.rotate_z(90)
        refmesh.translate((0, shift_y, 0))

        p = pv.Plotter(off_screen=True)
        p.add_mesh(refmesh, scalars="U", scalar_bar_args=sargs, cmap="coolwarm", clim=(vel_low,vel_high))
        p.add_mesh(resultmesh, scalars="U", scalar_bar_args=sargs, cmap="coolwarm", clim=(vel_low,vel_high))
        p.show(screenshot=output_U, cpos=(0, 0, 1), window_size=[4800, 4800])

    if casename == "reference":
        p = pv.Plotter(off_screen=True)
        p.add_mesh(resultmesh, scalars="p", scalar_bar_args=sargs, cmap="coolwarm", )
        p.show(screenshot=output_p, cpos=(0, 0, 1), window_size=[4800, 4800])
    else:
        shift_y = resultmesh.bounds[3] - resultmesh.bounds[2]
        refmesh = load_mesh(refpath)
        refmesh.rotate_z(90)
        refmesh.translate((0, shift_y, 0))

        p = pv.Plotter(off_screen=True)
        p.add_mesh(refmesh, scalars="p", scalar_bar_args=sargs, cmap="coolwarm", )
        p.add_mesh(resultmesh, scalars="p", scalar_bar_args=sargs, cmap="coolwarm", )
        p.show(screenshot=output_p, cpos=(0, 0, 1), window_size=[4800, 4800])

    if casename == "reference":
        p = pv.Plotter(off_screen=True)
        p.add_mesh(resultmesh, scalars="T", scalar_bar_args=sargs, cmap="coolwarm", clim=(temperature_low,temperature_high))
        p.show(screenshot=output_T, cpos=(0, 0, 1), window_size=[4800, 4800])
    else:
        shift_y = resultmesh.bounds[3] - resultmesh.bounds[2]
        refmesh = load_mesh(refpath)
        refmesh.rotate_z(90)
        refmesh.translate((0, shift_y, 0))

        p = pv.Plotter(off_screen=True)
        p.add_mesh(refmesh, scalars="T", scalar_bar_args=sargs, cmap="coolwarm", clim=(temperature_low,temperature_high))
        p.add_mesh(resultmesh, scalars="T", scalar_bar_args=sargs, cmap="coolwarm", clim=(temperature_low,temperature_high))
        p.show(screenshot=output_T, cpos=(0, 0, 1), window_size=[4800, 4800])

    if casename == "reference":
        p = pv.Plotter(off_screen=True)
        p.add_mesh(resultmesh, scalars="rho", scalar_bar_args=sargs, cmap="coolwarm", clim=(rho_low,rho_high))
        p.show(screenshot=output_rho, cpos=(0, 0, 1), window_size=[4800, 4800])
    else:
        shift_y = resultmesh.bounds[3] - resultmesh.bounds[2]
        refmesh = load_mesh(refpath)
        refmesh.rotate_z(90)
        refmesh.translate((0, shift_y, 0))

        p = pv.Plotter(off_screen=True)
        p.add_mesh(refmesh, scalars="rho", scalar_bar_args=sargs, cmap="coolwarm", clim=(rho_low,rho_high))
        p.add_mesh(resultmesh, scalars="rho", scalar_bar_args=sargs, cmap="coolwarm", clim=(rho_low,rho_high))
        p.show(screenshot=output_rho, cpos=(0, 0, 1), window_size=[4800, 4800])


    if casename == "reference":
        p = pv.Plotter(off_screen=True)
        p.add_mesh(resultmesh, scalars="TurbulentEnergyKinetic", scalar_bar_args=sargs, cmap="coolwarm", clim=(tke_low,tke_high))
        p.show(screenshot=output_tke, cpos=(0, 0, 1), window_size=[4800, 4800])
    else:
        shift_y = resultmesh.bounds[3] - resultmesh.bounds[2]
        refmesh = load_mesh(refpath)
        refmesh.rotate_z(90)
        refmesh.translate((0, shift_y, 0))

        p = pv.Plotter(off_screen=True)
        p.add_mesh(refmesh, scalars="TurbulentEnergyKinetic", scalar_bar_args=sargs, cmap="coolwarm", clim=(tke_low,tke_high))
        p.add_mesh(resultmesh, scalars="TurbulentEnergyKinetic", scalar_bar_args=sargs, cmap="coolwarm", clim=(tke_low,tke_high))
        p.show(screenshot=output_tke, cpos=(0, 0, 1), window_size=[4800, 4800])

def check_trace_residuals(input, output):
    def readCSV(path, start_row, header, seperator):
        df = pd.read_csv(path, skiprows=start_row, header=header, sep=seperator)
        return df

    converged = {}
    for i in input:
        residual = readCSV(i, 4, None, seperator="\s+")
        rmsMax = 1e-3
        if not all([True if i > rmsMax else False for i in residual[3][-100:]]):
            converged[input] = True
            print(i + " converged")
        else:
            converged[input] = False

            print(i + " not converged")
    if all(converged.values()):
        f = open(output[0], "a")
        f.write("OK")
        f.close()


def boundary_layer_values(input, output, casename, turbprod, relative_height, chordlength, ymin, ymax):
    mesh = load_mesh(input)
    mesh = mesh.extract_geometry()
    mesh = mesh.compute_normals()
    bounds = mesh.bounds

    zspan = np.sign((bounds[5] - bounds[4])) * (bounds[5] - bounds[4]) * float(relative_height) + bounds[4]
    slice = mesh.slice(normal="z", origin=(0, 0, zspan))

    xmin = min(slice.points[::, 0])
    xmax = max(slice.points[::, 0])
    #    slice.points -= np.array([xmin,0,0])
    displacement_thickness = slice["ThicknessDisplacement"]
    momentum_thickness = slice["ThicknessMomentum"]
    h12 = displacement_thickness / momentum_thickness
    xs = (slice.points[:, 0] - xmin) / xmax
    ys = h12

    plt.figure()
    scattersize = 2
    plt.scatter(xs, ys, s=scattersize, label=input)
    if casename != "reference":
        ref = input.replace(casename, "reference").replace(turbprod, "10")
        mesh = load_mesh(ref)
        mesh = mesh.extract_geometry()
        mesh = mesh.compute_normals()
        bounds = mesh.bounds
        slice = mesh.slice(normal="z", origin=(0, 0, zspan))
        xmin = min(slice.points[::, 0])
        xmax = max(slice.points[::, 0])
        #    slice.points -= np.array([xmin,0,0])
        displacement_thickness = slice["ThicknessDisplacement"]
        momentum_thickness = slice["ThicknessMomentum"]
        h12 = displacement_thickness / momentum_thickness
        xs = (slice.points[:, 0] - xmin) / xmax
        ys = h12
        plt.scatter(xs, ys, s=scattersize, label=ref, color="orange")

    plt.ylim((ymin, ymax))
    plt.xlabel("$c_{ax}$")
    plt.ylabel("$H_{12}$")
    plt.title(input)
    plt.grid()
    plt.legend()
    plt.savefig(output)
    plt.close()


def plot_wake_profile(input, output, xpos1, xpos2, casename, tprod):
    def sort_ylike(ys, us):
        return [list(t) for t in zip(*sorted(zip(ys, us)))]

    fig, axs = plt.subplots(1, 2, figsize=(13, 8))

    slice = load_mesh(input)
    line1 = slice.slice(normal="x", origin=(xpos1, 0, 0))

    ys_1 = line1.points[::, 1]
    us_1 = np.array([vecAbs(i) for i in line1["U"]])
    ys_1, us_1 = sort_ylike(ys_1, us_1)

    line2 = slice.slice(normal="x", origin=(xpos2, 0, 0))
    ys_2 = line2.points[::, 1]
    us_2 = np.array([vecAbs(i) for i in line2["U"]])
    ys_2, us_2 = sort_ylike(ys_2, us_2)

    axs[0].plot(ys_1, us_1, label=casename)
    axs[1].plot(ys_2, us_2, label=casename)

    if casename != "reference":
        ref = input.replace(casename, "reference").replace(tprod, "10")
        refslice = load_mesh(ref)
        rline1 = refslice.slice(normal="x", origin=(xpos1, 0, 0))
        rys_1 = rline1.points[::, 1]
        rus_1 = np.array([vecAbs(i) for i in rline1["U"]])
        rys_1, rus_1 = sort_ylike(rys_1, rus_1)

        rline2 = refslice.slice(normal="x", origin=(xpos2, 0, 0))
        rys_2 = rline2.points[::, 1]
        rus_2 = np.array([vecAbs(i) for i in rline2["U"]])
        rys_2, rus_2 = sort_ylike(rys_2, rus_2)

        axs[0].plot(rys_1, rus_1, label="reference")
        axs[1].plot(rys_2, rus_2, label="reference")

    fig.suptitle(input)
    axs[0].set_title("xpos1: " + str(xpos1))
    axs[1].set_title("xpos2: " + str(xpos2))
    axs[0].grid()
    axs[1].grid()
    plt.legend()
    # plt.title(input)
    plt.tight_layout()
    plt.savefig(output)
    plt.close()


def turbintensity_along_domain(input, output):
    mesh = load_mesh(input)

    res = 100
    bounds = mesh.bounds
    xstart = bounds[0]
    xend = bounds[1]
    xpos = (np.linspace(0, 1, res) * (xend - xstart) + xstart)[1:-1]

    TuLine = {"Tu": [], "x": []}
    for pos in xpos:
        xslice = mesh.slice(origin=(pos, 0, 0), normal="x")
        xslice = xslice.compute_normals()
        Tus = (xslice["TurbulentEnergyKinetic"] * 2 / 3) ** (.5) / np.array([vecAbs(i) for i in xslice["U"]])
        xslice["Tus"] = Tus
        xslice = xslice.point_data_to_cell_data()
        TuLine["Tu"].append(massflowAvePlane(xslice, "Tus"))
        TuLine["x"].append(pos)
    write_pickle(output, TuLine)


def turbintensity_along_domain_plots(input, output, name, tprod, lowlim, highlim):
    plt.figure()
    caseline = read_pickle(input)
    plt.plot(caseline["x"], caseline["Tu"], label=input)

    if name != "reference":
        ref = input.replace(name, "reference").replace(tprod, "10")
        refline = read_pickle(ref)
        plt.plot(refline["x"], refline["Tu"], label=ref)

    plt.ylim((lowlim, highlim))
    plt.legend()
    plt.grid()
    plt.title(input)
    plt.tight_layout()
    plt.savefig(output)
    plt.close()


def tke_along_domain_plots(input, output, name, tprod, xpos, lowlim, highlim):
    pv.set_plot_theme("document")

    res = 4800
    title_size = int(0.02 * res)
    sargs = dict(
        title_font_size=title_size,
        label_font_size=int(0.016 * res),
        shadow=True,
        n_labels=3,
        italic=True,
        # fmt="%.1f",
        font_family="arial",
    )

    p = pv.Plotter(off_screen=True)

    mesh = load_mesh(input)
    slice = mesh.slice(origin=(xpos, 0, 0), normal="x")
    slice.set_active_scalars("TurbulentEnergyKinetic")
    slice.rotate_z(90)
    p.add_mesh(slice, scalar_bar_args=sargs, cmap="coolwarm", clim=[lowlim, highlim])

    if name != "reference":
        ref = input.replace(name, "reference").replace(tprod, "10")
        refmesh = load_mesh(ref)
        refslice = refmesh.slice(origin=(xpos, 0, 0), normal="x")
        refslice.set_active_scalars("TurbulentEnergyKinetic")
        refslice.rotate_z(90)
        slicebounds = slice.bounds
        pitch = slicebounds[3] - slicebounds[2]
        refslice.translate((0, pitch, 0))
        p.add_mesh(refslice, scalar_bar_args=sargs, cmap="coolwarm", clim=[lowlim, highlim])

    #    p.set_focus(slice.center)
    #    p.camera_set = True
    # p.camera_position = [(-0.5428039203228628,-0.24600000000000039, -0.24600000000000039),
    #                     (0.12300000000000019, -0.3260885553590027, -0.12700000000000003),
    #                     (0.0, 0.0, 1.0)]
    # p.camera_set=True
    p.show(screenshot=output, window_size=[4800, 4800], title=input)

    # cpos missing


def plot_vol_differences(input_mesh, input_blade, output, name, tprod, cpos_param=None):
    mesh = load_mesh(input_mesh)
    blade = load_mesh(input_blade)
    pv.set_plot_theme("document")

    annotations = {
        0: "Unchanged",
    }

    res = 4800
    title_size = int(0.02 * res)
    sargs = dict(
        title_font_size=title_size,
        label_font_size=int(0.016 * res),
        shadow=True,
        n_labels=3,
        italic=True,
        # fmt="%.1f",
        font_family="arial",
    )
    differences = mesh.copy()
    if name != "reference":
        differences.clear_data()

        ref = input_mesh.replace(name, "reference").replace(str(tprod), "10")
        refmesh = load_mesh(ref)

        for an in ["TurbulentEnergyKinetic"]:
            differences[an] = -(refmesh[an] - mesh[an]) / refmesh[an]

        low = min(differences["TurbulentEnergyKinetic"])
        high = max(differences["TurbulentEnergyKinetic"])
        differences["TKE_change"] = differences["TurbulentEnergyKinetic"] / max([abs(high), abs(low)])
        maxdiff_id = np.argmax(np.abs(differences.point_data["TKE_change"]))
        maxdiff_point = differences.points[maxdiff_id]
        my_cmap_RdBu = create_transparent_cmap(pl.cm.seismic, low, high)

        differencespts = pv.PolyData(differences.points)
        differencespts["TKE_change"] = differences.point_data["TKE_change"]
        # edges = differences.extract_all_edges()
        feature_edges = mesh.extract_feature_edges()
        blade = blade.extract_geometry()
        blade.clear_data()
        p = pv.Plotter(off_screen=True)
        p.add_mesh(differencespts, annotations=annotations, scalars="TKE_change", cmap=my_cmap_RdBu,
                   scalar_bar_args=sargs, clim=[-1, 1], n_colors=1024)
        p.add_mesh(feature_edges, show_scalar_bar=False, color="black")
        p.add_mesh(blade, show_edges=True, opacity=0.1, color="black")
        p.add_mesh(maxdiff_point, color="black", point_size=40, render_points_as_spheres=True)
        if cpos_param:
            p.show(screenshot=output, window_size=[res, res], title=input_mesh, cpos=cpos_param)
        else:
            p.show(screenshot=output, window_size=[res, res], title=input_mesh)
    else:
        feature_edges = mesh.extract_feature_edges()
        p = pv.Plotter(off_screen=True)
        p.add_mesh(mesh, annotations=annotations, scalars="TurbulentEnergyKinetic", cmap=pl.cm.seismic,
                   scalar_bar_args=sargs)
        p.add_mesh(blade, show_edges=True, opacity=0.1)
        p.add_mesh(feature_edges, show_scalar_bar=False, color="black")

        if cpos_param:
            p.show(screenshot=output, window_size=[res, res], title=input_mesh, cpos=cpos_param)
        else:
            p.show(screenshot=output, window_size=[res, res], title=input_mesh)


def create_transparent_cmap(cmap, low, high):
    """
    :param cmap: matplotlib.pylab.cm - colormap
    :return: cap - transparent
    """
    # norm = max([abs(low),abs(high)])
    my_cmap = cmap(np.arange(cmap.N))
    # Set alpha
    my_cmap[:, -1] = abs(np.linspace(-1, 1, cmap.N))
    # Create new colormap
    my_cmap = ListedColormap(my_cmap)
    return my_cmap


def outflow_onedimensional(input, output, direction, facenormal):
    mesh = load_mesh(input)
    # todo: this is a workaround
    mesh.rotate_x(90)
    bounds = mesh.cell_centers().bounds
    dirdir = {"x": 0, "y": 2, "z": 4}
    d = dirdir[facenormal]
    interd = dirdir[direction]
    xlow = (d + 2) % 6
    xhigh = (d + 3) % 6
    ylow = (d + 4) % 6
    yhigh = (d + 5) % 6

    res = 1000
    stepx = (bounds[xhigh] - bounds[xlow]) / res
    stepy = (bounds[yhigh] - bounds[ylow]) / res
    xi = np.arange(bounds[xlow], bounds[xhigh] + stepx, stepx)
    yi = np.arange(bounds[ylow], bounds[yhigh] + stepy, stepy)
    xi, yi = np.meshgrid(xi, yi)

    xs, ys = mesh.points[::, int((d / 2 + 1) % 3)], mesh.points[::, int((d / 2 + 2) % 3)]

    interpolated_2d = {}
    for arrname in mesh.array_names:
        if len(mesh[arrname].shape) == 1:
            # mask = (xi > 0.5) & (xi < 0.6) & (yi > 0.5) & (yi < 0.6)
            zi = griddata((xs, ys), mesh[arrname], (xi, yi), method='linear')
        else:
            zi = np.array([griddata((xs, ys), mesh[arrname][::, i], (xi, yi), method='linear') for i in
                           range(len(mesh[arrname][0]))])
        interpolated_2d[arrname] = zi

    interpolated_1d = {}
    rhofield = interpolated_2d["rho"]

    for arrname in interpolated_2d.keys():
        array = interpolated_2d[arrname]
        dimension = array.ndim
        arr = np.zeros(len(array))
        if dimension == 3:
            arr = []
            for i in range(dimension):
                arrmean = np.average(array[i], axis=0, weights=rhofield)
                arr.append(arrmean)
            arr = np.dstack(arr)[0]
        else:
            arr = np.average(array, axis=0, weights=rhofield)
            #arr += arrmean
        interpolated_1d[arrname] = arr
    interpolated_1d["xs"] = xi
    interpolated_1d["ys"] = yi
    write_pickle(output, interpolated_1d)


def plot_outflow_onedimensional(input, outTKE, outU, outAlpha, name, prod_raw):
    onedim_dict = read_pickle(input)
    refpath = os.path.join(input.replace(name, "reference").replace(prod_raw, "10"))

    reference_dict = read_pickle(refpath)

    velocities = onedim_dict["U"]

    xs = onedim_dict["xs"]
    velocities_normed = np.array([vecAbs(i) for i in velocities])
    alphas = np.array([radToGrad(vecAngle(i, np.array([1, 0, 0]))) for i in velocities])

    tke = onedim_dict["TurbulentEnergyKinetic"]

    if name != "reference":

        velocities_ref = reference_dict["U"]
        velocities_ref_normed = np.array([vecAbs(i) for i in velocities_ref])
        alphas_ref_normed = np.array([radToGrad(vecAngle(i, np.array([1, 0, 0]))) for i in velocities_ref])
        tke_ref = reference_dict["TurbulentEnergyKinetic"]

        plt.figure()
        plt.plot(velocities_normed, color="blue", label=input)
        plt.plot(velocities_ref_normed, color="red", label=refpath)
        plt.grid(True, linestyle='-', linewidth=1)
        plt.legend()
        plt.title(input)
        plt.savefig(outU)
        plt.close()

        plt.figure()
        plt.plot(tke, color="blue")
        plt.plot(tke_ref, color="red")
        plt.grid(True, linestyle='-', linewidth=1)
        plt.legend()
        plt.title(input)
        plt.savefig(outTKE)
        plt.close()

        plt.figure()
        plt.plot(alphas, color="blue", label=input)
        plt.plot(alphas_ref_normed, color="red", label=refpath)
        plt.grid(True, linestyle='-', linewidth=1)
        plt.legend()
        plt.title(input)
        plt.savefig(outAlpha)
        plt.close()


    else:
        plt.figure()
        plt.plot(velocities_normed, color="blue", label=input)
        plt.grid(True, linestyle='-', linewidth=1)
        plt.legend()
        plt.title(input)
        plt.savefig(outU)
        plt.close()

        plt.figure()
        plt.plot(tke, color="blue")
        plt.grid(True, linestyle='-', linewidth=1)
        plt.legend()
        plt.title(input)
        plt.savefig(outTKE)
        plt.close()

        plt.figure()
        plt.plot(alphas, color="blue", label=input)
        plt.grid(True, linestyle='-', linewidth=1)
        plt.legend()
        plt.title(input)
        plt.savefig(outAlpha)
        plt.close()
