import numpy as np
import pyvista as pv

from NTR.utils.geom_functions import extract_geo_paras, calcMidPoints, calcMidPassageStreamLine, rotate_points
from NTR.utils.externals.tecplot.tecplot_functions import writeTecplot1DFile


def create_geometry(path_profile_coords, x_inlet, x_outlet, pitch, unit, blade_shift, alpha,midline_tol, verbose=False):
    # =============================================================================
    # Daten Einlesen
    # =============================================================================
    points = np.loadtxt(path_profile_coords)
    unitcoeff = 0
    if unit == "m":
        unitcoeff = 1
    elif unit == "mm":
        unitcoeff = 1/1000
    points *= unitcoeff



    # =============================================================================
    # Bestimmung Profilparameter
    # =============================================================================
    #x_ss, y_ss, x_ps, y_ps = sortProfilePoints(x_raw, y_raw, alpha=alpha)
    psPoly,ssPoly,ind_vk,ind_hk, midsPoly, beta_meta_01, beta_meta_02 = extract_geo_paras(points,alpha,midline_tol)


    #x_mids, y_mids = calcMidPoints(x_ss, y_ss, x_ps, y_ps)
    #x_mids, y_mids = zip(*sorted(zip(x_mids, y_mids)))
    x_mids = midsPoly.points[::,0]
    y_mids = midsPoly.points[::,1]
    x_ss = ssPoly.points[::,0]
    y_ss = ssPoly.points[::,1]
    x_ps = psPoly.points[::,0]
    y_ps = psPoly.points[::,1]

    stagger_angle = 0#np.rad2deg(np.arctan((y_mids[-1] - y_mids[-0]) / (x_mids[-1] - x_mids[-0])))
    #x_ss, y_ss = rotate_points([0, 0], x_ss, y_ss, - stagger_angle)
    #x_ps, y_ps = rotate_points([0, 0], x_ps, y_ps, -stagger_angle)
    #x_mids, y_mids = calcMidPoints(x_ss, y_ss, x_ps, y_ps)
    #x_mids, y_mids = zip(*sorted(zip(x_mids, y_mids)))
    x_mpsl, y_mpsl = calcMidPassageStreamLine(x_mids, y_mids, (beta_meta_01 - 90) - stagger_angle - 90,
                                              (beta_meta_02 - 90) - stagger_angle - 90, x_inlet, x_outlet,
                                              pitch * np.sin(np.deg2rad(beta_meta_01)))
    #x_mpsl, y_mpsl = rotate_points([0, 0], x_mpsl, y_mpsl, + stagger_angle)
    #x_ss, y_ss = rotate_points([0, 0], x_ss, y_ss, +stagger_angle)
    #x_ps, y_ps = rotate_points([0, 0], x_ps, y_ps, +stagger_angle)
    #x_mids, y_mids = rotate_points([0, 0], x_mids, y_mids, +stagger_angle)
    y_upper = np.array(y_mpsl) + blade_shift # +0.05*pitch
    y_lower = y_upper - pitch
    writeTecplot1DFile('01_Meshing/geom.dat', ['x', 'z'], ['druckseite', 'saugseite', 'lower peri', 'upper peri', 'skelett'],
                       [[x_ss, y_ss], [x_ps, y_ps], [x_mpsl, y_lower], [x_mpsl, y_upper], [x_mids[::-1], y_mids[::-1]]],
                       'obere Kurvenverlauf des Kanals')

    if verbose:
        plotter = pv.Plotter()
        psPoly = pv.PolyData(np.stack((x_ss, y_ss, np.zeros(len(x_ss)))).T)
        ssPoly = pv.PolyData(np.stack((x_ps, y_ps, np.zeros(len(x_ps)))).T)
        mpslPolyLow = pv.PolyData(np.stack((x_mpsl, y_lower, np.zeros(len(x_mpsl)))).T)
        mpslPolyUpper = pv.PolyData(np.stack((x_mpsl, y_upper, np.zeros(len(x_mpsl)))).T)
        midsPoly = pv.PolyData(np.stack((x_mids, y_mids, np.zeros(len(x_mids)))).T)

        plotter.add_mesh(psPoly)
        plotter.add_mesh(ssPoly)
        plotter.add_mesh(mpslPolyLow)
        plotter.add_mesh(mpslPolyUpper)
        plotter.add_mesh(midsPoly)
        plotter.show_axes()
        plotter.show()
