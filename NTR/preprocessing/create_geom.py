import numpy as np
import pyvista as pv
import os

from NTR.utils.geom_functions.pointcloud import calcConcaveHull
from NTR.utils.geom_functions.profileparas import calcMidPassageStreamLine, extract_vk_hk, extractSidePolys, \
    midline_from_sides, angles_from_mids
from NTR.utils.externals.tecplot.tecplot_functions import writeTecplot1DFile
from NTR.utils.filehandling import write_pickle, yaml_dict_read
from NTR.utils.geom_functions.pyvista_utils import lines_from_points


def create_geometry(path_profile_coords, x_inlet, x_outlet, pitch, unit, blade_shift, alpha, span_z, casepath,
                    verbose=False):
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
                                                                                                           verbose)

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

    if True:
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

        plotter.add_mesh(psPoly,color="red")
        plotter.add_mesh(ssPoly,color="blue")
        plotter.add_mesh(mpslPolyLow)
        plotter.add_mesh(mpslPolyUpper)
        plotter.add_mesh(midsPoly,point_size=5)

        plotter.add_mesh(psPoly_asline)
        plotter.add_mesh(ssPoly_asline)
        plotter.add_mesh(mpslPolyLow_asline)
        plotter.add_mesh(mpslPolyUpper_asline)
        plotter.add_mesh(midsPoly_asline)

        plotter.show_axes()
        plotter.show()


def extract_geo_paras(points, alpha, verbose):
    """
    This function is extracting profile-data as stagger-angle, midline, psPoly, ssPoly and more from a set of points
    Be careful, you need a suitable alpha-parameter in order to get the right geometry
    The calculation of the leading-edge and trailing-edge index needs time and its not 100% reliable (yet)
    Keep in mind, to check the results!
    :param points: array of points in 3d with the shape (n,3)
    :param alpha: nondimensional alpha-coefficient (calcConcaveHull)
    :param verbose: bool for plots
    :return: points, psPoly, ssPoly, ind_vk, ind_hk, midsPoly, camber_angle_vk, camber_angle_hk
    """

    origPoly = pv.PolyData(points)
    xs, ys = calcConcaveHull(points[:, 0], points[:, 1], alpha)
    points = np.stack((xs, ys, np.zeros(len(xs)))).T
    sortedPoly = pv.PolyData(points)

    ind_hk, ind_vk, veronoi_mid = extract_vk_hk(origPoly, sortedPoly)
    psPoly, ssPoly = extractSidePolys(ind_hk, ind_vk, sortedPoly)
    midsPoly = midline_from_sides(ind_hk, ind_vk, points, psPoly, ssPoly)
    camber_angle_hk, camber_angle_vk = angles_from_mids(midsPoly)

    if verbose:
        p = pv.Plotter()
        p.add_mesh(points, color="orange", label="points")
        p.add_mesh(psPoly, color="green", label="psPoly")
        p.add_mesh(ssPoly, color="black", label="ssPoly")
        p.add_mesh(midsPoly, color="black", label="midsPoly")
        p.add_mesh(veronoi_mid, color="yellow", label="veronoi_mid")
        p.add_legend()
        p.show()

    return points, psPoly, ssPoly, ind_vk, ind_hk, midsPoly, camber_angle_vk, camber_angle_hk


def run_create_geometry_frompointcloud(settings_yaml):
    case_path = os.path.abspath(os.path.dirname(settings_yaml))
    settings = yaml_dict_read(settings_yaml)
    meshpath = os.path.join(case_path, "01_Meshing")
    if not os.path.isdir(meshpath):
        os.mkdir(meshpath)
    print(os.path.abspath(case_path))

    print("create_geometry")
    ptstxtfile = os.path.join(os.path.abspath(case_path), settings["geom"]["ptcloud_profile"])
    outpath = os.path.join(os.path.dirname(os.path.abspath(settings_yaml)),"00_Ressources","01_Geometry")
    create_geometry(ptstxtfile,
                    settings["geom"]["x_inlet"],
                    settings["geom"]["x_outlet"],
                    settings["geometry"]["pitch"],
                    settings["geom"]["ptcloud_profile_unit"],
                    settings["geom"]["shift_domain"],
                    settings["geometry"]["alpha"],
                    settings["mesh"]["extrudeLength"],
                    outpath,)
    return 0
