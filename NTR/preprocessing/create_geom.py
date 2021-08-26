import numpy as np
import pyvista as pv
import os

from NTR.utils.geom_functions.pointcloud import calcConcaveHull
from NTR.utils.geom_functions.profileparas import calcMidPassageStreamLine, extract_vk_hk, extractSidePolys, \
    midline_from_sides, angles_from_mids
from NTR.utils.externals.tecplot.tecplot_functions import writeTecplot1DFile
from NTR.utils.filehandling import write_pickle, yaml_dict_read
from NTR.utils.geom_functions.pyvista_utils import lines_from_points
from NTR.database.naca_airfoil_creator import naca
from NTR.utils.mathfunctions import vecAbs


def create_geometry_frompointcloud(path_profile_coords, x_inlet, x_outlet, pitch, unit, blade_shift, alpha, span_z,
                                   casepath,
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

        plotter.add_mesh(psPoly, color="red")
        plotter.add_mesh(ssPoly, color="blue")
        plotter.add_mesh(mpslPolyLow)
        plotter.add_mesh(mpslPolyUpper)
        plotter.add_mesh(midsPoly, point_size=5)

        plotter.add_mesh(psPoly_asline)
        plotter.add_mesh(ssPoly_asline)
        plotter.add_mesh(mpslPolyLow_asline)
        plotter.add_mesh(mpslPolyUpper_asline)
        plotter.add_mesh(midsPoly_asline)

        plotter.show_axes()
        plotter.show()
    return geo_dict



def create_geometry_fromnacaairfoil(nacadigits, numberofpoints, finite_TE, half_cosine_spacing, x_inlet, x_outlet,
                                    pitch, camberlength, blade_shift, staggerangle, span_z, casepath, verbose=False):
    # =============================================================================
    # Bestimmung Profilparameter
    # =============================================================================
    sortedPoints, psPoly, ssPoly, ind_vk, ind_hk, midsPoly, beta_meta_01, beta_meta_02 = create_naca_geoparas(
        nacadigits,
        numberofpoints,
        finite_TE,
        half_cosine_spacing,
        camberlength,
        staggerangle,
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

    geo_dict = {"points": sortedPoints,
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

        plotter.add_mesh(psPoly, color="red")
        plotter.add_mesh(ssPoly, color="blue")
        plotter.add_mesh(mpslPolyLow)
        plotter.add_mesh(mpslPolyUpper)
        plotter.add_mesh(midsPoly, point_size=5)

        plotter.add_mesh(psPoly_asline)
        plotter.add_mesh(ssPoly_asline)
        plotter.add_mesh(mpslPolyLow_asline)
        plotter.add_mesh(mpslPolyUpper_asline)
        plotter.add_mesh(midsPoly_asline)

        plotter.show_axes()
        plotter.show()
    return geo_dict


def create_naca_geoparas(nacadigits, numberofpoints, finite_TE, half_cosine_spacing, camberlength, staggerangle, verbose=False):
    ptsx, ptsy = naca(nacadigits, numberofpoints,finite_TE, half_cosine_spacing)
    ind_hk = 0
    ind_vk = numberofpoints
    points = np.stack((ptsx[:-1], ptsy[:-1], np.zeros(numberofpoints * 2) - 1)).T
    scalegeo = camberlength / vecAbs(points[ind_hk] - points[ind_vk])
    poly = pv.PolyData(points)
    poly.points*=scalegeo
    poly.rotate_z(staggerangle)
    points = poly.points
    ssPoly, psPoly = extractSidePolys(ind_hk, ind_vk, poly)
    midsPoly = midline_from_sides(ind_hk, ind_vk, points, psPoly, ssPoly)
    camber_angle_hk, camber_angle_vk = angles_from_mids(midsPoly)

    if verbose:
        p = pv.Plotter()
        p.add_mesh(points, color="orange", label="points")
        p.add_mesh(psPoly, color="green", label="psPoly")
        p.add_mesh(ssPoly, color="black", label="ssPoly")
        p.add_mesh(midsPoly, color="black", label="midsPoly")
        p.add_legend()
        p.show()
    return points, psPoly, ssPoly, ind_vk, ind_hk, midsPoly, camber_angle_vk, camber_angle_hk


def extract_geo_paras(points, alpha, verbose=False):
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
    print("starting extract_geo_paras")
    print()
    print("sorting point-data...")

    origPoly = pv.PolyData(points)
    xs, ys = calcConcaveHull(points[:, 0], points[:, 1], alpha)
    points = np.stack((xs, ys, np.zeros(len(xs)))).T
    sortedPoly = pv.PolyData(points)

    print("extracting LE/TE-indices...")
    ind_hk, ind_vk, veronoi_mid = extract_vk_hk(origPoly, sortedPoly, verbose)
    print("extracting ps/ss-sides...")
    psPoly, ssPoly = extractSidePolys(ind_hk, ind_vk, sortedPoly)
    print("calculating midline...")
    midsPoly = midline_from_sides(ind_hk, ind_vk, points, psPoly, ssPoly)
    print("calculating angles...")
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


def run_create_geometry(settings_yaml):
    case_path = os.path.abspath(os.path.dirname(settings_yaml))
    settings = yaml_dict_read(settings_yaml)
    meshpath = os.path.join(case_path, "01_Meshing")
    if not os.path.isdir(meshpath):
        os.mkdir(meshpath)
    print(os.path.abspath(case_path))

    print("create_geometry")
    if settings["geom"]["algorithm"] == "from_pointcloud":
        ptstxtfile = os.path.join(os.path.abspath(case_path), settings["geom"]["ptcloud_profile"])
        outpath = os.path.join(os.path.dirname(os.path.abspath(settings_yaml)), "00_Ressources", "01_Geometry")
        geo_dict = create_geometry_frompointcloud(ptstxtfile,
                                       settings["geom"]["x_inlet"],
                                       settings["geom"]["x_outlet"],
                                       settings["geometry"]["pitch"],
                                       settings["geom"]["ptcloud_profile_unit"],
                                       settings["geom"]["shift_domain"],
                                       settings["geometry"]["alpha"],
                                       settings["mesh"]["extrudeLength"],
                                       outpath, )

    if settings["geom"]["algorithm"] == "naca_airfoil_generator":
        outpath = os.path.join(os.path.dirname(os.path.abspath(settings_yaml)), "00_Ressources", "01_Geometry")
        geo_dict = create_geometry_fromnacaairfoil(settings["geom"]["naca_digits"],
                                        settings["geom"]["numberofpoints"],
                                        settings["geom"]["finite_TE"],
                                        settings["geom"]["half_cosine_spacing"],
                                        settings["geom"]["x_inlet"],
                                        settings["geom"]["x_outlet"],
                                        settings["geometry"]["pitch"],
                                        settings["geom"]["camberlength"],
                                        settings["geom"]["shift_domain"],
                                        settings["geom"]["staggerangle"],
                                        settings["mesh"]["extrudeLength"],
                                        outpath, )
    #blockpoints(geo_dict,settings)
    return 0

"""
def blockpoints(geodict,settings):
    ssPoly = geodict["sidePolys"]["ssPoly"]
    psPoly = geodict["sidePolys"]["psPoly"]
    ylower = geodict["periodics"]["ylower"]


    return 0
"""
