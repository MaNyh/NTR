# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 20:33:27 2019

@author: Mark Ziesse / Malte Nyhuis
"""


import matplotlib.pyplot as plt
import os
import numpy as np
import pyvista as pv

from NTR.utils.geom_functions import sortProfilePoints, calcMidPassageStreamLine, calcMeanCamberLine, getBoundaryValues, getGeom2DVTUSLice2


def createProbesProfileDict(case, pden_Probes_Profile_SS, pden_Probes_Profile_PS,
                            interval_time_steps_probes, output_path, tolerance=1e-6):

    """
    :param path_blade_surface: vtk-mesh
    :param pden_Probes_Profile_SS: integer (slicer)
    :param pden_Probes_Profile_PS: integer (slicer)
    :param interval_time_steps_probes: integer
    :param output_path: pathstring
    :param tolerance: float
    :return: openFoamDict
    """
    # blade_surface einlesen und Normalenvektor bestimmen der
    blade_surface = case.mesh_loaded_dict["blade"]
    blade_surface = blade_surface.compute_normals()

    bladebounds = blade_surface.bounds
    midspan_z = (bladebounds[5]-bladebounds[4])/2

    cut_plane = blade_surface.slice(normal="z", origin=(0, 0, midspan_z))

    # Mittelschnitt erstellen

    points = cut_plane.points
    surface_normals = cut_plane.cell_normals
    # Punkte extrahieren
    x_values = []
    y_values = []

    for i in range(len(points)):
        point = points[i]
        normal = surface_normals[i]

        x_values.append(point[0] - tolerance * normal[0])
        y_values.append(point[1] - tolerance * normal[1])


    # Nach Durck und Saugseite sortieren
    x_ss, y_ss, x_ps, y_ps = sortProfilePoints(x_values, y_values, case.CascadeCoeffs.alpha)

    x_bl_ss = x_ss[::int(pden_Probes_Profile_SS)]
    y_bl_ss = y_ss[::int(pden_Probes_Profile_SS)]

    x_bl_ps = x_ps[::int(pden_Probes_Profile_PS)]
    y_bl_ps = y_ps[::int(pden_Probes_Profile_PS)]

    z_bl_ss = []
    z_bl_ps = []

    for i in range(len(x_bl_ss)):
        z_bl_ss.append(midspan_z)

    for i in range(len(x_bl_ps)):
        z_bl_ps.append(midspan_z)

    data_file = open(os.path.join(output_path, 'Probes_Profile_Dict.txt'), 'w')

    data_file.write("""    Probes_Profile
    {
        type                probes;
        libs                ("libsampling.so");
        writeControl        timeStep;
        writeInterval       """ + str(int(interval_time_steps_probes)) + """;

            fields
            (
                U
                p
                T
                gradU
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


def createProbesStreamlineDict(case, nop_Probes_Streamline, save_dir,
                               interval_time_steps_probes, beta_01, beta_02, teilung):

    x_bounds, y_bounds, x_profil, y_profil , midspan_z = getGeom2DVTUSLice2(case)

    y_inlet, x_inlet, y_outlet, x_outlet, x_lower_peri, y_lower_peri, x_upper_peri, y_upper_peri = getBoundaryValues(
        x_bounds, y_bounds)

    x_mids, y_mids, x_ss, y_ss, x_ps, y_ps, x_vk, y_vk, x_hk, y_hk = calcMeanCamberLine(x_profil, y_profil, beta_01,
                                                                                        beta_02)

    x_mpsl, y_mpsl = calcMidPassageStreamLine(x_mids, y_mids, beta_01, beta_02, max(x_inlet), min(x_outlet), teilung)

    #x_probes = []
    #y_probes = []
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

    data_file = open(os.path.join(save_dir, 'Probes_Streamline_Dict.txt'), 'w')

    data_file.write("""
Probes_Streamline
{
type                probes;
libs                ("libsampling.so");
writeControl        timeStep;
writeInterval       """ + str(int(interval_time_steps_probes)) + """;

fields
(
    U
    p
    T
    gradU
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
    plt.plot(x_lower_peri, y_lower_peri, '-y', lw=1, label='lower_peri')
    plt.plot(x_upper_peri, y_upper_peri, '-c', lw=1, label='upper_peri')
    plt.plot(x_probes, y_probes, 'x', label='Probes_Streamline', color="darkorange")
    plt.legend(loc='best')
    plt.savefig(os.path.join(save_dir, 'kontrollplot_probes_streamline.pdf'))
    plt.close('all')



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


