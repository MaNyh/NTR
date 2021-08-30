# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 20:33:27 2019

@author: Mark Ziesse / Malte Nyhuis
"""

import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
import os
import numpy as np
import pyvista as pv

from NTR.utils.geom_functions.geom_utils import getBoundaryValues, equi_points
from NTR.utils.geom_functions.spline import refine_spline
from NTR.utils.geom_functions.pyvista_utils import polyline_from_points

def createProbesProfileDict(geoparas_dict, midspan_z, pden_Probes_Profile_SS, pden_Probes_Profile_PS,
                            interval_time_steps_probes, output_path,  start_time, end_time, tolerance):
    """
    :param path_blade_surface: vtk-mesh
    :param pden_Probes_Profile_SS: integer (slicer)
    :param pden_Probes_Profile_PS: integer (slicer)
    :param interval_time_steps_probes: integer
    :param output_path: pathstring
    :param tolerance: float
    :return: openFoamDict
    """
    ssPoly = geoparas_dict["sidePolys"]["ssPoly"]
    ref_ss_x, ref_ss_y = refine_spline(ssPoly.points[::, 0], ssPoly.points[::, 1], 4000)
    ref_ss_points = np.stack((ref_ss_x, ref_ss_y, np.zeros(len(ref_ss_y)))).T
    ref_ssPoly = pv.PolyData(ref_ss_points)
    ref_ss_poly = polyline_from_points(ref_ssPoly.points)
    ref_ss_face = ref_ss_poly.extrude((0, 0, midspan_z * 2)).compute_normals()
    ref_ss_face_shift = ref_ss_face.copy()
    ref_ss_face_shift.points += tolerance * ref_ss_face_shift.point_arrays["Normals"]
    ref_ss_cut = ref_ss_face_shift.slice(normal="z", origin=(0, 0, midspan_z))

    psPoly = geoparas_dict["sidePolys"]["psPoly"]
    ref_ps_x, ref_ps_y = refine_spline(psPoly.points[::, 0], psPoly.points[::, 1], 4000)
    ref_ps_points = np.stack((ref_ps_x, ref_ps_y, np.zeros(len(ref_ps_y)))).T
    ref_psPoly = pv.PolyData(ref_ps_points)
    ref_ps_poly = polyline_from_points(ref_psPoly.points)
    ref_ps_face = ref_ps_poly.extrude((0, 0, midspan_z * 2)).compute_normals()
    ref_ps_face_shift = ref_ps_face.copy()
    ref_ps_face_shift.points += tolerance * ref_ps_face_shift.point_arrays["Normals"]
    ref_ps_cut = ref_ps_face_shift.slice(normal="z", origin=(0, 0, midspan_z))

    x_ss_shift = ref_ss_cut.points[::,0]
    y_ss_shift = ref_ss_cut.points[::,1]
    z_bl_ss = ref_ps_cut.points[::,2]

    x_ps_shift = ref_ps_cut.points[::,0]
    y_ps_shift = ref_ps_cut.points[::,1]
    z_bl_ps = ref_ps_cut.points[::,2]

    x_bl_ss, y_bl_ss = refine_spline(x_ss_shift, y_ss_shift, pden_Probes_Profile_SS)
    x_bl_ps, y_bl_ps = refine_spline(x_ps_shift, y_ps_shift, pden_Probes_Profile_PS)



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
    """
    plt.close('all')
    plt.figure(figsize=(8, 8))

    plt.plot(x_bl_ss, y_bl_ss, 'xr', lw=1, label='probes ss_profile')
    plt.plot(x_bl_ps, y_bl_ps, 'xb', lw=1, label='probes ps_profile')

    plt.legend(loc='best')
    plt.savefig(os.path.join(os.path.abspath(output_path), 'kontrollplot_probes_profile.pdf'))
    plt.close('all')
    """
    outprobes = {"probes pressure-side": np.stack((x_bl_ss, y_bl_ss)),
                 "probes suction-side": np.stack((x_bl_ps, y_bl_ps)),
                 }

    return outprobes


def createProbesStreamlineDict(nop_Probes_Streamline, save_dir,
                               interval_time_steps_probes, start_time, end_time, geoparas_dict):


    midspan_z = geoparas_dict["span_z"]/2
    y_lower = geoparas_dict["periodics"]["ylower"]
    y_upper = geoparas_dict["periodics"]["yupper"]

    outer = y_lower.merge(y_upper)

    x_bounds = outer.points[::,0]
    y_bounds = outer.points[::,1]

    y_inlet, x_inlet, y_outlet, x_outlet, x_lower_peri, y_lower_peri, x_upper_peri, y_upper_peri = getBoundaryValues(
        x_bounds, y_bounds)

    x_mpsl = np.array(geoparas_dict["midpassagestreamLine"]["x_mpsl"])
    y_mpsl = np.array(geoparas_dict["midpassagestreamLine"]["y_mpsl"])

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
    """
    plt.close('all')
    plt.figure(figsize=(8, 8))
    plt.plot(x_inlet, y_inlet, '-r', lw=1, label='inlet')
    plt.plot(x_outlet, y_outlet, '-b', lw=1, label='outlet')
    plt.plot(*geoparas_dict["points"][::,0:2], '.k', lw=1, label='profil')
    plt.plot(x_lower_peri, y_lower_peri, '-y', lw=1, label='lower_peri')
    plt.plot(x_upper_peri, y_upper_peri, '-c', lw=1, label='upper_peri')
    plt.plot(x_probes, y_probes, 'x', label='Probes_Streamline', color="darkorange")
    plt.legend(loc='best')
    plt.savefig(os.path.join(save_dir, 'kontrollplot_probes_streamline.pdf'))
    plt.close('all')
    """
    outprobes = {"probes streamline": np.stack((x_probes, y_probes)),
                 }

    return outprobes


def createProbesInletOutlet(geodat_dict, interval_time_steps_probes, output_path, start_time, end_time):

    midspan_z = geodat_dict["span_z"]/2
    y_lower = geodat_dict["periodics"]["ylower"]
    y_upper = geodat_dict["periodics"]["yupper"]

    outer = y_lower.merge(y_upper)

    x_bounds = outer.points[::,0]
    y_bounds = outer.points[::,1]

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


def createXSliceProbes(geodat_dict, nop, x_slice_1, x_slice_2, interval_time_steps_probes,
                       output_path, start_time, end_time):

    midspan_z = geodat_dict["span_z"]/2
    y_lower = geodat_dict["periodics"]["ylower"]
    y_upper = geodat_dict["periodics"]["yupper"]

    perbounds = y_lower.merge(y_upper)

    slice_1 = perbounds.slice(origin=(x_slice_1, 0, 0), normal=(1, 0, 0))
    slice_2 = perbounds.slice(origin=(x_slice_2, 0, 0), normal=(1, 0, 0))

    ys_1 = slice_1.points[::, 1]
    ys_2 = slice_2.points[::, 1]

    y1max = max(ys_1)
    y1min = min(ys_1)

    y2max = max(ys_2)
    y2min = min(ys_2)

    dy_shift = (y2max - y2min) / nop / 2
    y1_probes = np.linspace(y1min, y1max, nop, endpoint=False) + dy_shift
    y2_probes = np.linspace(y2min, y2max, nop, endpoint=False) + dy_shift

    data_file = open(os.path.join(output_path, 'Probes_XSlices_Dict'), 'w')

    data_file.write("""    Probes_XSlices
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


def create_vk_stagflow_probes(geom_paras, nop, length, angle, interval_time_steps_probes,
                              output_path, start_time, end_time):

    vk_point = geom_paras["sortedPoly"][geom_paras["hk_vk_idx"]["ind_vk"]]
    stagnationLine = pv.Line((0, 0, 0), (-length, 0, 0), nop - 1)
    stagnationLine.rotate_z(angle)
    stagnationLine.translate(vk_point)
    midspan_z = geom_paras["span_z"] / 2
    x_probes = stagnationLine.points[::, 0]
    y_probes = stagnationLine.points[::, 1]
    z_probes = stagnationLine.points[::, 2] + midspan_z
    data_file = open(os.path.join(output_path, 'Probes_VKstagnation_Dict'), 'w')

    data_file.write("""
Probes_VKstagnation
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

    outprobes = {"Probes_VK_StagnationLine": np.stack((x_probes, y_probes)),
                 }

    return outprobes


def create_of_les_probe_dicts(case_settings, geo_ressources):

    midspan_z = case_settings["mesh"]["extrudeLength"] / 2

    output_path = os.path.join("02_Simcase", "system")

    probes = {}

    if case_settings["probing"]["probes"]["profile_probing"]:
        sampling_rate = case_settings["probes"]["profile_probes"]["sampling_rate"]
        timestepinterval = int(float(sampling_rate) ** -1 / float(case_settings["openfoam_cascade_les_settings"]["timestep"]))
        outprobes = createProbesProfileDict(geo_ressources,
                                            midspan_z,
                                            case_settings["probes"]["profile_probes"]["pden_ps"],
                                            case_settings["probes"]["profile_probes"]["pden_ss"],
                                            timestepinterval,
                                            output_path,
                                            case_settings["probes"]["profile_probes"]["start_time"],
                                            case_settings["probes"]["profile_probes"]["end_time"],
                                            case_settings["probes"]["profile_probes"]["tolerance"]
                                            )
        for k, v in outprobes.items():
            probes[k] = v

    if case_settings["probing"]["probes"]["streamline_probing"]:

        sampling_rate = case_settings["probes"]["streamline_probes"]["sampling_rate"]
        timestepinterval = int(float(sampling_rate) ** -1 / float(case_settings["openfoam_cascade_les_settings"]["timestep"]))
        outprobes = createProbesStreamlineDict(case_settings["probes"]["streamline_probes"]["nop_streamline"],
                                               output_path,
                                               timestepinterval,
                                               case_settings["probes"]["streamline_probes"]["start_time"],
                                               case_settings["probes"]["streamline_probes"]["end_time"],
                                               geo_ressources)

        for k, v in outprobes.items():
            probes[k] = v

    if case_settings["probing"]["probes"]["inletoutletvelocity_probing"]:
        sampling_rate = case_settings["probes"]["inletoutlet_probes"]["sampling_rate"]
        timestepinterval = int(float(sampling_rate) ** -1 / float(case_settings["openfoam_cascade_les_settings"]["timestep"]))
        outprobes = createProbesInletOutlet(geo_ressources,
                                            timestepinterval,
                                            output_path,
                                            case_settings["probes"]["inletoutlet_probes"]["start_time"],
                                            case_settings["probes"]["inletoutlet_probes"]["end_time"])
        for k, v in outprobes.items():
            probes[k] = v

    if case_settings["probing"]["probes"]["xslice_probing"]:

        sampling_rate = case_settings["probes"]["xsclicing_probes"]["sampling_rate"]
        timestepinterval = int(float(sampling_rate) ** -1 / float(case_settings["openfoam_cascade_les_settings"]["timestep"]))
        outprobes = createXSliceProbes(geo_ressources,
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
        timestepinterval = int(float(sampling_rate) ** -1 / float(case_settings["openfoam_cascade_les_settings"]["timestep"]))
        create_inletoutletave_probe_dict(case_settings["probes"]["inletoutletfieldave_probing"]["start_time"],
                                         case_settings["probes"]["inletoutletfieldave_probing"]["end_time"],
                                         timestepinterval,
                                         output_path)

    if case_settings["probing"]["probes"]["vk_stagnationflow_probing"]:
        sampling_rate = case_settings["probes"]["inletoutletfieldave_probing"]["sampling_rate"]

        U = [float(i) for i in case_settings["case"]["case_parameters"]["U"]["UINLET"].split(" ")]
        angle = np.arccos(U[0] / U[1]) * 180 / np.pi
        timestepinterval = int(float(sampling_rate) ** -1 / float(case_settings["openfoam_cascade_les_settings"]["timestep"]))

        outprobes = create_vk_stagflow_probes(geo_ressources,
                                              case_settings["probes"]["vk_stagnationflow_probing"]["nop"],
                                              case_settings["probes"]["vk_stagnationflow_probing"]["length"],
                                              angle,
                                              timestepinterval,
                                              output_path,
                                              case_settings["probes"]["inletoutletfieldave_probing"]["start_time"],
                                              case_settings["probes"]["inletoutletfieldave_probing"]["end_time"]
                                              )
        for k, v in outprobes.items():
            probes[k] = v

    #midspan_z = geo_ressources["span_z"]/2
    y_lower = geo_ressources["periodics"]["ylower"]
    y_upper = geo_ressources["periodics"]["yupper"]
    ssPoly = geo_ressources["sidePolys"]["ssPoly"]
    psPoly = geo_ressources["sidePolys"]["psPoly"]
    inlet = geo_ressources["flowbounds"]["inletPoly"]
    outlet = geo_ressources["flowbounds"]["outletPoly"]

    plt.close('all')
    plt.figure(figsize=(8, 8))
    plt.plot(inlet.points[::,0],inlet.points[::,1], '-r', lw=1, label='inlet')
    plt.plot(outlet.points[::,0],outlet.points[::,1], '-b', lw=1, label='outlet')
    plt.plot(ssPoly.points[::,0],ssPoly.points[::,1], '-k', lw=1, label='ss_profil')
    plt.plot(psPoly.points[::,0],psPoly.points[::,1], '-k', lw=1, label='ps_profil')
    plt.plot(y_lower.points[::,0],y_lower.points[::,1], '-y', lw=1, label='lower_peri')
    plt.plot(y_upper.points[::,0],y_upper.points[::,1], '-c', lw=1, label='upper_peri')

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
