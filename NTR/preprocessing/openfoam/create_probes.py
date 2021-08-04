# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 20:33:27 2019

@author: Mark Ziesse / Malte Nyhuis
"""

import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
import os
import numpy as np

from NTR.utils.geom_functions import sortProfilePoints, calcMidPassageStreamLine, calcMeanCamberLine, getBoundaryValues, \
    getGeom2DVTUSLice2
from NTR.utils.pyvista_utils import load_mesh


def createProbesProfileDict(blade_surface, pden_Probes_Profile_SS, pden_Probes_Profile_PS,
                            interval_time_steps_probes, output_path, alpha, start_time, end_time, tolerance):
    """
    :param path_blade_surface: vtk-mesh
    :param pden_Probes_Profile_SS: integer (slicer)
    :param pden_Probes_Profile_PS: integer (slicer)
    :param interval_time_steps_probes: integer
    :param output_path: pathstring
    :param tolerance: float
    :return: openFoamDict
    """

    bladebounds = blade_surface.bounds
    blade_surface = blade_surface.compute_normals()
    surface_normals = blade_surface.point_arrays["Normals"]

    midspan_z = (bladebounds[5] - bladebounds[4]) / 2

    cut_plane = blade_surface.slice(normal="z", origin=(0, 0, midspan_z))

    # Mittelschnitt erstellen

    points = cut_plane.points
    phelp = points[:, [0, 1]]

    #
    # Punkte extrahieren

    # Nach Durck und Saugseite sortieren
    x_ss, y_ss, x_ps, y_ps = sortProfilePoints(points[::, 0], points[::, 1], alpha)

    x_ss_shift = []
    y_ss_shift = []

    x_ps_shift = []
    y_ps_shift = []

    for pt in np.stack([x_ss, y_ss]).T:
        idx = np.where(pt == phelp)[0][0]
        point = points[idx]
        normal = surface_normals[idx]

        x_ss_shift.append(point[0] + tolerance * normal[0])
        y_ss_shift.append(point[1] + tolerance * normal[1])

    for pt in np.stack([x_ps, y_ps]).T:
        idx = np.where(pt == phelp)[0][0]
        point = points[idx]
        normal = surface_normals[idx]

        x_ps_shift.append(point[0] - tolerance * normal[0])
        y_ps_shift.append(point[1] - tolerance * normal[1])

    x_bl_ss = x_ss_shift[::int(pden_Probes_Profile_SS)]
    y_bl_ss = y_ss_shift[::int(pden_Probes_Profile_SS)]

    x_bl_ps = x_ps_shift[::int(pden_Probes_Profile_PS)]
    y_bl_ps = y_ps_shift[::int(pden_Probes_Profile_PS)]

    z_bl_ss = []
    z_bl_ps = []

    for i in range(len(x_bl_ss)):
        z_bl_ss.append(midspan_z)

    for i in range(len(x_bl_ps)):
        z_bl_ps.append(midspan_z)

    data_file = open(os.path.join(output_path, 'Probes_Profile_Dict'), 'w')

    data_file.write("""    Probes_Profile
    {
        type                probes;
        libs                ("libsampling.so");
        writeControl        timeStep;
        writeInterval       """ + str(int(interval_time_steps_probes)) + """;
        timeStart           """ + str(start_time) + """;
        timeEnd             """ + str(end_time) + """;

        fields
            (
                U
                p
                T
                nut
                rho
            );

        // number of probes: """ + str(len(x_bl_ss) + len(x_bl_ps)) + """

        probeLocations
        (\n""")

    data_file.write('\t\t\t//Probes auf der Saugseite ' + str(len(x_bl_ss)) + '\n\n')

    for i in range(len(x_bl_ss)):
        data_file.write('\t\t\t(' + str(x_bl_ss[i]) + '\t\t' + str(y_bl_ss[i]) + '\t\t' + str(z_bl_ss[i]) + ')\n')

    data_file.write('\n\t\t\t//Probes auf der Druckseite ' + str(len(x_bl_ps)) + '\n\n')

    for i in range(len(x_bl_ps)):
        data_file.write('\t\t\t(' + str(x_bl_ps[i]) + '\t\t' + str(y_bl_ps[i]) + '\t\t' + str(z_bl_ps[i]) + ')\n')

    data_file.write("""        );
\t}""")
    data_file.close()

    plt.close('all')
    plt.figure(figsize=(8, 8))

    plt.plot(x_bl_ss, y_bl_ss, 'xr', lw=1, label='probes ss_profile')
    plt.plot(x_bl_ps, y_bl_ps, 'xb', lw=1, label='probes ss_profile')

    plt.legend(loc='best')
    plt.savefig(os.path.join(os.path.abspath(output_path), 'kontrollplot_probes_profile.pdf'))
    plt.close('all')

    outprobes = {"probes pressure-side": np.stack((x_bl_ss, y_bl_ss)),
                 "probes suction-side": np.stack((x_bl_ps, y_bl_ps)),
                 }

    return outprobes


def createProbesStreamlineDict(mesh, alpha, nop_Probes_Streamline, save_dir,
                               interval_time_steps_probes,start_time, end_time, geoparas_dict):

    x_bounds, y_bounds, x_profil, y_profil, midspan_z = getGeom2DVTUSLice2(mesh, alpha)

    y_inlet, x_inlet, y_outlet, x_outlet, x_lower_peri, y_lower_peri, x_upper_peri, y_upper_peri = getBoundaryValues(
        x_bounds, y_bounds)

    x_mpsl = np.array(geoparas_dict["midpassagestreamLine"]).T[::,0]
    y_mpsl = np.array(geoparas_dict["midpassagestreamLine"]).T[::,1]

    # x_probes = []
    # y_probes = []
    z_probes = []

    nop = int(nop_Probes_Streamline)

    xn, yn = equi_points(x_mpsl, y_mpsl, nop)
    for i in range(nop):
        z_probes.append(midspan_z)

    x_probes = xn
    y_probes = yn

    dist = np.sqrt((x_probes[0] - x_probes[1]) ** 2 + (y_probes[0] - y_probes[1]) ** 2)

    x_probes[0] = x_probes[0] + 0.00001 * dist
    x_probes[-1] = x_probes[-1] - 0.00001 * dist

    data_file = open(os.path.join(save_dir, 'Probes_Streamline_Dict'), 'w')

    data_file.write("""
Probes_Streamline
{
type                probes;
libs                ("libsampling.so");
writeControl        timeStep;
writeInterval       """ + str(int(interval_time_steps_probes)) + """;
timeStart           """ + str(start_time) + """;
timeEnd             """ + str(end_time) + """;

fields
(
    U
    p
    T
    nut
    rho
);

// number of probes: """ + str(nop) + """

probeLocations
(\n""")

    for i in range(len(x_probes)):
        data_file.write('\t(' + str(x_probes[i]) + '\t\t' + str(y_probes[i]) + '\t\t' + str(z_probes[i]) + ')\n')

    data_file.write("""        );
}""")
    data_file.close()

    plt.close('all')
    plt.figure(figsize=(8, 8))
    plt.plot(x_inlet, y_inlet, '-r', lw=1, label='inlet')
    plt.plot(x_outlet, y_outlet, '-b', lw=1, label='outlet')
    plt.plot(x_profil, y_profil, '.k', lw=1, label='profil')
    plt.plot(x_lower_peri, y_lower_peri, '-y', lw=1, label='lower_peri')
    plt.plot(x_upper_peri, y_upper_peri, '-c', lw=1, label='upper_peri')
    plt.plot(x_probes, y_probes, 'x', label='Probes_Streamline', color="darkorange")
    plt.legend(loc='best')
    plt.savefig(os.path.join(save_dir, 'kontrollplot_probes_streamline.pdf'))
    plt.close('all')

    outprobes = {"probes streamline": np.stack((x_probes, y_probes)),
                 }

    return outprobes


def equi_points(x, y, nop):
    M = 10000

    x_new = np.linspace(min(x), max(x), M)
    y_new = np.interp(x_new, x, y)

    # berechnet die laenge der Stromlinie

    l_sl = 0

    for i in range(len(x_new)):
        if i > 0:
            l_sl = l_sl + np.sqrt((x_new[i] - x_new[i - 1]) ** 2 + (y_new[i] - y_new[i - 1]) ** 2)

    xn = []
    yn = []

    dist = l_sl / (nop - 1)

    l_p = 0

    for i in range(len(x_new)):
        if i > 0:
            l_p = l_p + np.sqrt((x_new[i] - x_new[i - 1]) ** 2 + (y_new[i] - y_new[i - 1]) ** 2)

            if l_p >= dist and i != nop - 1:
                xn.append(x_new[i])
                yn.append(y_new[i])
                l_p = 0

        if i == 0:
            xn.append(x_new[i])
            yn.append(y_new[i])

        if i == len(x_new) - 1:
            xn.append(x_new[-1])
            yn.append(y_new[-1])

    return xn, yn


def createProbesInletOutlet(mesh, alpha, interval_time_steps_probes, output_path, start_time, end_time):
    x_bounds, y_bounds, x_profil, y_profil, midspan_z = getGeom2DVTUSLice2(mesh, alpha)

    y_inlet, x_inlet, y_outlet, x_outlet, x_lower_peri, y_lower_peri, x_upper_peri, y_upper_peri = getBoundaryValues(
        x_bounds, y_bounds)

    x_probe_inlet = (max(x_inlet) + min(x_inlet)) / 2
    y_probe_inlet = (max(y_inlet) + min(y_inlet)) / 2

    x_probe_outlet = (max(x_outlet) + min(x_outlet)) / 2
    y_probe_outlet = (max(y_outlet) + min(y_outlet)) / 2

    data_file = open(os.path.join(output_path, 'Probes_InletOutlet_Dict'), 'w')

    data_file.write("""    Probes_Profile
        {
            type                probes;
            libs                ("libsampling.so");
            writeControl        timeStep;
            writeInterval       """ + str(int(interval_time_steps_probes)) + """;
            timeStart           """ + str(start_time) + """;
            timeEnd             """ + str(end_time) + """;

                fields
                (
                    U
                    p
                );

            probeLocations
            (\n""")

    data_file.write('\t\t\t//Probe Inlet ' + '\n\n')
    data_file.write('\t\t\t(' + str(x_probe_inlet) + '\t\t' + str(y_probe_inlet) + '\t\t' + str(midspan_z) + ')\n')
    data_file.write('\t\t\t//Probe Outlet ' + '\n\n')
    data_file.write('\t\t\t(' + str(x_probe_outlet) + '\t\t' + str(y_probe_outlet) + '\t\t' + str(midspan_z) + ')\n')

    data_file.write("""        );
    \t}""")
    data_file.close()

    outprobes = {"probes inlet": np.stack((x_probe_inlet, y_probe_inlet)),
                 "probes outlet": np.stack((x_probe_outlet, y_probe_outlet)),
                 }

    return outprobes


def createXSliceProbes(mesh, nop, x_slice_1, x_slice_2, interval_time_steps_probes, output_path, start_time, end_time):
    bounds = mesh.bounds
    midspan_z = (bounds[5] + bounds[4]) / 2
    slice_1 = mesh.slice(origin=(x_slice_1, 0, midspan_z), normal=(1, 0, 0))
    line_1 = slice_1.slice(normal=(0, 0, 1))
    slice_2 = mesh.slice(origin=(x_slice_2, 0, midspan_z), normal=(1, 0, 0))
    line_2 = slice_2.slice(normal=(0, 0, 1))

    ys_1 = line_1.points[::, 1]
    ys_2 = line_2.points[::, 1]

    y1max = max(ys_1)
    y1min = min(ys_1)

    y2max = max(ys_2)
    y2min = min(ys_2)

    y1_probes = np.linspace(y1min, y1max, nop, endpoint=True)
    y2_probes = np.linspace(y2min, y2max, nop, endpoint=True)

    data_file = open(os.path.join(output_path, 'Probes_XSlices_Dict'), 'w')

    data_file.write("""    Probes_Profile
        {
            type                probes;
            libs                ("libsampling.so");
            writeControl        timeStep;
            writeInterval       """ + str(int(interval_time_steps_probes)) + """;
            timeStart           """ + str(start_time) + """;
            timeEnd             """ + str(end_time) + """;

                fields
                (
                    U
                );

            // number of probes: """ + str(len(y1_probes) + len(y2_probes)) + """

            probeLocations
            (\n""")

    data_file.write('\t\t\t//Probes auf X-Slice 1 ' + str(len(y1_probes)) + '\n\n')

    for i in range(len(y1_probes)):
        data_file.write('\t\t\t(' + str(x_slice_1) + '\t\t' + str(y1_probes[i]) + '\t\t' + str(midspan_z) + ')\n')

    data_file.write('\n\t\t\t//Probes auf X-Slice 2 ' + str(len(y2_probes)) + '\n\n')

    for i in range(len(y2_probes)):
        data_file.write('\t\t\t(' + str(x_slice_2) + '\t\t' + str(y2_probes[i]) + '\t\t' + str(midspan_z) + ')\n')

    data_file.write("""        );
    \t}""")
    data_file.close()

    outprobes = {"probes xslice1": np.stack((x_slice_1 * np.ones(len(y1_probes)), y1_probes)),
                 "probes xslice2": np.stack((x_slice_2 * np.ones(len(y1_probes)), y2_probes)),
                 }

    return outprobes


def create_probe_dicts(case_settings,geo_ressources):

    domain = load_mesh(case_settings["probing"]["domain"])
    blade = load_mesh(case_settings["probing"]["blade"])
    alpha = case_settings["geometry"]["alpha"]
    beta_01 = geo_ressources["beta_metas"][0]
    beta_02 = geo_ressources["beta_metas"][1]
    pitch = case_settings["geometry"]["pitch"]

    output_path = case_settings["probing"]["output_path"]

    probes = {}

    if case_settings["probing"]["probes"]["profile_probing"]:
        sampling_rate = case_settings["probes"]["profile_probes"]["sampling_rate"]
        timestepinterval = int(float(sampling_rate) ** -1 / float(case_settings["case_settings"]["timestep"]))
        outprobes = createProbesProfileDict(blade,
                                            case_settings["probes"]["profile_probes"]["pden_ps"],
                                            case_settings["probes"]["profile_probes"]["pden_ss"],
                                            timestepinterval,
                                            output_path,
                                            alpha,
                                            case_settings["probes"]["profile_probes"]["start_time"],
                                            case_settings["probes"]["profile_probes"]["end_time"],
                                            case_settings["probes"]["profile_probes"]["tolerance"]
                                            )
        for k, v in outprobes.items():
            probes[k] = v

    if case_settings["probing"]["probes"]["streamline_probing"]:

        sampling_rate = case_settings["probes"]["streamline_probes"]["sampling_rate"]
        timestepinterval = int(float(sampling_rate) ** -1 / float(case_settings["case_settings"]["timestep"]))
        outprobes = createProbesStreamlineDict(domain,
                                               alpha,
                                               case_settings["probes"]["streamline_probes"]["nop_streamline"],
                                               output_path,
                                               timestepinterval,
                                               case_settings["probes"]["streamline_probes"]["start_time"],
                                               case_settings["probes"]["streamline_probes"]["end_time"],
                                               geo_ressources)

        for k, v in outprobes.items():
            probes[k] = v

    if case_settings["probing"]["probes"]["inletoutletvelocity_probing"]:
        sampling_rate = case_settings["probes"]["inletoutlet_probes"]["sampling_rate"]
        timestepinterval = int(float(sampling_rate) ** -1 / float(case_settings["case_settings"]["timestep"]))
        outprobes = createProbesInletOutlet(domain, alpha,
                                            timestepinterval,
                                            output_path,
                                            case_settings["probes"]["inletoutlet_probes"]["start_time"],
                                            case_settings["probes"]["inletoutlet_probes"]["end_time"])
        for k, v in outprobes.items():
            probes[k] = v

    if case_settings["probing"]["probes"]["xslice_probing"]:

        sampling_rate = case_settings["probes"]["xsclicing_probes"]["sampling_rate"]
        timestepinterval = int(float(sampling_rate) ** -1 / float(case_settings["case_settings"]["timestep"]))
        outprobes = createXSliceProbes(domain,
                                       case_settings["probes"]["xsclicing_probes"]["nop"],
                                       case_settings["probes"]["xsclicing_probes"]["x_slice_one"],
                                       case_settings["probes"]["xsclicing_probes"]["x_slice_two"],
                                       timestepinterval,
                                       output_path,
                                       case_settings["probes"]["xsclicing_probes"]["start_time"],
                                       case_settings["probes"]["xsclicing_probes"]["end_time"])
        for k, v in outprobes.items():
            probes[k] = v

    if case_settings["probing"]["probes"]["inletoutletfieldave_probing"]:
        sampling_rate = case_settings["probes"]["inletoutletfieldave_probing"]["sampling_rate"]
        timestepinterval = int(float(sampling_rate) ** -1 / float(case_settings["case_settings"]["timestep"]))
        create_inletoutletave_probe_dict(case_settings["probes"]["inletoutletfieldave_probing"]["start_time"],
                                         case_settings["probes"]["inletoutletfieldave_probing"]["end_time"],
                                         timestepinterval,
                                         output_path)

    x_bounds, y_bounds, x_profil, y_profil, midspan_z = getGeom2DVTUSLice2(domain, alpha)

    y_inlet, x_inlet, y_outlet, x_outlet, x_lower_peri, y_lower_peri, x_upper_peri, y_upper_peri = getBoundaryValues(
        x_bounds, y_bounds)

    plt.close('all')
    plt.figure(figsize=(8, 8))
    plt.plot(x_inlet, y_inlet, '-r', lw=1, label='inlet')
    plt.plot(x_outlet, y_outlet, '-b', lw=1, label='outlet')
    plt.plot(x_profil, y_profil, '.k', lw=1, label='profil')
    plt.plot(x_lower_peri, y_lower_peri, '-y', lw=1, label='lower_peri')
    plt.plot(x_upper_peri, y_upper_peri, '-c', lw=1, label='upper_peri')

    allowed_colors = list(dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS).keys())
    del allowed_colors[allowed_colors.index("white")]
    del allowed_colors[allowed_colors.index("whitesmoke")]
    del allowed_colors[allowed_colors.index("seashell")]
    del allowed_colors[allowed_colors.index("linen")]
    del allowed_colors[allowed_colors.index("floralwhite")]
    del allowed_colors[allowed_colors.index("azure")]
    del allowed_colors[allowed_colors.index("mintcream")]
    del allowed_colors[allowed_colors.index("ghostwhite")]
    del allowed_colors[allowed_colors.index("aliceblue")]
    del allowed_colors[allowed_colors.index("lavenderblush")]
    del allowed_colors[allowed_colors.index("lightcyan")]
    del allowed_colors[allowed_colors.index("oldlace")]
    del allowed_colors[allowed_colors.index("antiquewhite")]
    del allowed_colors[allowed_colors.index("mistyrose")]

    for k, v in probes.items():
        color = allowed_colors.pop(np.random.randint(len(allowed_colors) - 1))

        plt.plot(v[0], v[1], color, marker="X", lw=1, markersize=3, label=k)

    plt.legend(loc='best')
    plt.savefig(os.path.join('kontrollplot_probes.pdf'))
    plt.close('all')


def create_inletoutletave_probe_dict(start_time, end_time, interval_time_steps_probes, output_path):
    with open(os.path.join(output_path, 'Probes_inletoutletave_Dict'), 'w') as data_file:
        data_file.write("""

MassflowInlet
    {
        type                surfaceFieldValue;
        libs                ("libfieldFunctionObjects.so");
        writeControl        timeStep;
        writeInterval       """ + str(int(interval_time_steps_probes)) + """;
        timeStart           """ + str(start_time) + """;
        timeEnd             """ + str(end_time) + """;

        log                     true;
        writeFields             false;
        regionType              patch;
        name                    INLET;
        operation               sum;

            fields
            (
                phi
            );
    }


AverValuesInlet
    {
        type                    surfaceFieldValue;
        libs                    ("libfieldFunctionObjects.so");
        writeControl            timeStep;
        writeInterval           """ + str(int(interval_time_steps_probes)) + """;
        timeStart           """ + str(start_time) + """;
        timeEnd             """ + str(end_time) + """;
        log                     true;
        writeFields             false;
        regionType              patch;
        name                    INLET;
        operation               areaAverage;

            fields
            (
                U
                p
                rho
                T

            );
    }


MassflowOutlet
    {
        type                surfaceFieldValue;
        libs                ("libfieldFunctionObjects.so");
        writeControl        timeStep;
        writeInterval       """ + str(int(interval_time_steps_probes)) + """;
        timeStart           """ + str(start_time) + """;
        timeEnd             """ + str(end_time) + """;

        log                     true;
        writeFields             false;
        regionType              patch;
        name                    OUTLET;
        operation               sum;

            fields
            (
                phi
            );
    }


AverValuesOutlet
    {
        type                    surfaceFieldValue;
        libs                    ("libfieldFunctionObjects.so");
        writeControl            timeStep;
        writeInterval       """ + str(int(interval_time_steps_probes)) + """;
        timeStart           """ + str(start_time) + """;
        timeEnd             """ + str(end_time) + """;
        log                     true;
        writeFields             false;
        regionType              patch;
        name                    OUTLET;
        operation               areaAverage;

            fields
            (
                U
                p
                rho
                T

            );
    }
    """)
    return 0
