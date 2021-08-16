import numpy as np
import pyvista as pv
import os

from NTR.utils.geom_functions import extract_geo_paras, calcMidPassageStreamLine
from NTR.utils.externals.tecplot.tecplot_functions import writeTecplot1DFile
from NTR.utils.filehandling import write_pickle
from NTR.utils.pyvista_utils import lines_from_points


def create_geometry(path_profile_coords, x_inlet, x_outlet, pitch, unit, blade_shift, alpha, midline_tol, span_z,
                    casepath,
                    verbose=True):
    # =============================================================================
    # Daten Einlesen
    # =============================================================================
    points = np.loadtxt(path_profile_coords)
    unitcoeff = 0
    if unit == "m":
        unitcoeff = 1
    elif unit == "mm":
        unitcoeff = 1 / 1000
    points *= unitcoeff

    # =============================================================================
    # Bestimmung Profilparameter
    # =============================================================================
    sortedPoints, psPoly, ssPoly, ind_vk, ind_hk, midsPoly, beta_meta_01, beta_meta_02 = extract_geo_paras(points,
                                                                                                           alpha,
                                                                                                           midline_tol)

    x_mids = midsPoly.points[::, 0]
    y_mids = midsPoly.points[::, 1]
    x_ss = ssPoly.points[::, 0]
    y_ss = ssPoly.points[::, 1]
    x_ps = psPoly.points[::, 0]
    y_ps = psPoly.points[::, 1]

    stagger_angle = np.rad2deg(np.arctan((y_mids[-1] - y_mids[-0]) / (x_mids[-1] - x_mids[-0])))

    x_mpsl, y_mpsl = calcMidPassageStreamLine(x_mids, y_mids, beta_meta_01, beta_meta_02,
                                              x_inlet, x_outlet, pitch)

    y_upper = np.array(y_mpsl) + blade_shift
    per_y_upper = pv.lines_from_points(np.stack((np.array(x_mpsl),
                                                 np.array(y_upper),
                                                 np.zeros(len(x_mpsl)))).T)
    y_lower = y_upper - pitch
    per_y_lower = pv.lines_from_points(np.stack((np.array(x_mpsl),
                                                 np.array(y_lower),
                                                 np.zeros(len(x_mpsl)))).T)

    inlet_pts = np.array([per_y_lower.points[per_y_lower.points[::, 0].argmin()],
                          per_y_upper.points[per_y_upper.points[::, 0].argmin()]])

    inletPoly = pv.Line(*inlet_pts)
    outlet_pts = np.array([per_y_lower.points[per_y_lower.points[::, 0].argmax()],
                           per_y_upper.points[per_y_upper.points[::, 0].argmax()]])

    outletPoly = pv.Line(*outlet_pts)

    writeTecplot1DFile('01_Meshing/geom.dat', ['x', 'z'],
                       ['druckseite', 'saugseite', 'lower peri', 'upper peri', 'skelett'],
                       [[x_ps, y_ps], [x_ss, y_ss], [x_mpsl, y_lower], [x_mpsl, y_upper], [x_mpsl[::-1], y_mpsl[::-1]]],
                       'obere Kurvenverlauf des Kanals')

    geo_dict = {"points": points,
                "sortedPoly": sortedPoints,
                "sidePolys": {"psPoly": psPoly, "ssPoly": ssPoly},
                "hk_vk_idx": {"ind_vk": ind_vk, "ind_hk": ind_hk},
                "periodics": {"ylower": per_y_lower, "yupper": per_y_upper},
                "flowbounds": {"inletPoly": inletPoly, "outletPoly": outletPoly},
                "midsPoly": midsPoly,
                "beta_metas": {"beta_meta_01": beta_meta_01, "beta_meta_02": beta_meta_02},
                "stagger_angle": stagger_angle,
                "midpassagestreamLine": {"x_mpsl": x_mpsl, "y_mpsl": y_mpsl},
                "xpos_in_out": {"x_inlet": x_inlet, "x_outlet": x_outlet},
                "pitch": pitch,
                "span_z": span_z
                }

    geo_filename = "geometry.pkl"
    write_pickle(os.path.join(casepath, geo_filename), geo_dict)

    if verbose:
        plotter = pv.Plotter()
        psPoly = pv.PolyData(np.stack((x_ss, y_ss, np.zeros(len(x_ss)))).T)
        ssPoly = pv.PolyData(np.stack((x_ps, y_ps, np.zeros(len(x_ps)))).T)
        mpslPolyLow = pv.PolyData(np.stack((x_mpsl, y_lower, np.zeros(len(x_mpsl)))).T)
        mpslPolyUpper = pv.PolyData(np.stack((x_mpsl, y_upper, np.zeros(len(x_mpsl)))).T)
        midsPoly = pv.PolyData(np.stack((x_mids, y_mids, np.zeros(len(x_mids)))).T)

        psPoly_asline = lines_from_points(np.stack((x_ss, y_ss, np.zeros(len(x_ss)))).T)
        ssPoly_asline = lines_from_points(np.stack((x_ps, y_ps, np.zeros(len(x_ps)))).T)
        mpslPolyLow_asline = lines_from_points(np.stack((x_mpsl, y_lower, np.zeros(len(x_mpsl)))).T)
        mpslPolyUpper_asline = lines_from_points(np.stack((x_mpsl, y_upper, np.zeros(len(x_mpsl)))).T)
        midsPoly_asline = lines_from_points(np.stack((x_mids, y_mids, np.zeros(len(x_mids)))).T)

        plotter.add_mesh(psPoly)
        plotter.add_mesh(ssPoly)
        plotter.add_mesh(mpslPolyLow)
        plotter.add_mesh(mpslPolyUpper)
        plotter.add_mesh(midsPoly)

        plotter.add_mesh(psPoly_asline)
        plotter.add_mesh(ssPoly_asline)
        plotter.add_mesh(mpslPolyLow_asline)
        plotter.add_mesh(mpslPolyUpper_asline)
        plotter.add_mesh(midsPoly_asline)

        plotter.show_axes()
        plotter.show()
