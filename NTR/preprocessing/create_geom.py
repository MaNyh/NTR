import numpy as np
import pyvista as pv
import os

from NTR.utils.geom_functions.pointcloud import calcConcaveHull
from NTR.utils.geom_functions.profileparas import calcMidPassageStreamLine, extract_vk_hk, extractSidePolys, \
    midline_from_sides, angles_from_mids
from NTR.utils.externals.tecplot.tecplot_functions import writeTecplot1DFile
from NTR.utils.filehandling import write_pickle, yaml_dict_read
from NTR.utils.pyvista_utils import lines_from_points, plot_geometry_tofile
from NTR.database.data_generators.naca_airfoil_creator import naca
from NTR.utils.mathfunctions import vecAbs
from NTR.preprocessing.prep import prep_geo
from NTR.database.case_dirstructure import casedirs


def create_geometry_frompointcloud(path_profile_coords, settings, casepath, verbose=False):
    x_inlet = settings["geometry"]["x_inlet"]
    x_outlet = settings["geometry"]["x_outlet"]
    pitch = settings["geometry"]["pitch"]
    unit = settings["geometry"]["ptcloud_profile_unit"]
    blade_shift = settings["geometry"]["shift_domain"]
    alpha = settings["geometry"]["alpha"]
    span_z = settings["mesh"]["extrudeLength"]

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
    sortedPoints, psPoly, ssPoly, ind_vk, ind_hk, midsPoly, beta_meta_01, beta_meta_02, camber_angle = extract_geo_paras(
        points,
        alpha,
        verbose)

    x_mids = midsPoly.points[::, 0]
    y_mids = midsPoly.points[::, 1]
    x_ss = ssPoly.points[::, 0]
    y_ss = ssPoly.points[::, 1]
    x_ps = psPoly.points[::, 0]
    y_ps = psPoly.points[::, 1]

    stagger_angle = np.rad2deg(np.arctan((y_mids[-1] - y_mids[-0]) / (x_mids[-1] - x_mids[-0])))
    chordlength = vecAbs(midsPoly.points[0] - midsPoly.points[-1])
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

    writeTecplot1DFile(os.path.join(casepath, '01_Meshing', 'geom.dat'), ['x', 'z'],
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
    """
    # this might be used later on for the definition of a generic blocking-algorithm
    noplines = 101
    ssLinePts_X, ssLinePts_Y = refine_spline(ssPoly.points[::, 0], ssPoly.points[::, 1], noplines)
    ssLine = lines_from_points(np.stack((ssLinePts_X, ssLinePts_Y, np.zeros(noplines))).T)

    psLinePts_X, psLinePts_Y = refine_spline(psPoly.points[::, 0], psPoly.points[::, 1], noplines)
    psLine = lines_from_points(np.stack((psLinePts_X, psLinePts_Y, np.zeros(noplines))).T)

    inlet_ylow = inlet_pts[0]  ## pt9
    inlet_southprofile = inlet_ylow + settings["mesh"]["yPerLowHGridBlockPitchStart"] * settings["geometry"][
        "pitch"] * np.array([0, 1, 0])  # pt11
    inlet_northprofile = inlet_ylow + settings["mesh"]["yPerHighHGridBlockPitchStart"] * settings["geometry"][
        "pitch"] * np.array([0, 1, 0])  # pt12

    vkplane_ylow = 1
    hkplane_ylow = 1
    outlet_ylow = outlet_pts[0]  # pt17
    outlet_southprofile = outlet_ylow + settings["mesh"]["yPerLowHGridBlockPitchStart"] * settings["geometry"][
        "pitch"] * np.array([0, 1, 0])  # pt18
    outlet_northprofile = outlet_ylow + settings["mesh"]["yPerHighHGridBlockPitchStart"] * settings["geometry"][
        "pitch"] * np.array([0, 1, 0])  # pt19
    midx_ysouth_2d = np.array(refine_spline(midsPoly.points[:, 0], midsPoly.points[:, 1], 3))[1]
    midx_ysouth_z = 0
    midx_ysouth = np.array([midx_ysouth_2d[0], midx_ysouth_2d[1], midx_ysouth_z])+ settings["mesh"]["yPerHighHGridBlockPitchStart"] * settings["geometry"][
        "pitch"] * np.array([0, 1, 0])
    inlet_yhigh = inlet_pts[1]  # pt13
    vkplane_yhigh = refine_spline(per_y_upper.points[:,0],per_y_upper.points[:,1],noplines)
    vkplane_yhigh = np.stack((vkplane_yhigh[0],vkplane_yhigh[1],np.zeros(noplines))).T
    vkplane_yhigh = vkplane_yhigh[int(100 * settings["mesh"]["yPerHighHGridBlockPitchStart"])]

    hkplane_yhigh = 1
    outlet_yhigh = outlet_pts[1]  # pt21

    inlet_profile_ylow = 1
    vkplane_profile_ylow = 1
    hkplane_profile_ylow = 1
    outlet_profile_ylow = 1

    inlet_profile_yhigh = 1
    vkplane_profile_yhigh = 1
    hkplane_profile_yhigh = 1
    outlet_profile_yhigh = 1

    vk_profileblockpt_yhigh = 1
    vk_profileblockpt_yhigh = 1

    vk_profileblockpt_ylow = 1
    vk_profileblockpt_ylow = 1

    psPoly_line = lines_from_points(psPoly.points)
    psPoly_surf = psPoly_line.extrude(vector=(0, 0, 0.1))
    psPoly_surf = psPoly_surf.compute_normals()
    psPoly_surf.points += 0.01 * psPoly_surf.point_arrays["Normals"]
    psPoly_surf.translate((0, 0, -0.05))
    psPoly_slice = psPoly_surf.slice(origin=(0, 0, 0), normal=(0, 0, 1))

    ssPoly_line = lines_from_points(ssPoly.points)
    ssPoly_surf = ssPoly_line.extrude(vector=(0, 0, 0.1))
    ssPoly_surf = ssPoly_surf.compute_normals()
    ssPoly_surf.points += 0.01 * ssPoly_surf.point_arrays["Normals"]
    ssPoly_surf.translate((0, 0, -0.05))
    ssPoly_slice = ssPoly_surf.slice(origin=(0, 0, 0), normal=(0, 0, 1))
    p = pv.Plotter()
    p.add_mesh(psPoly, label="psPoly", color="red", point_size=5)
    p.add_mesh(ssPoly, label="ssPoly", color="green", point_size=5)
    p.add_mesh(psPoly_slice)
    p.add_mesh(ssPoly_slice)
    p.add_legend()
    p.add_axes()
    p.show()
    """
    geo_filename = "geometry.pkl"
    write_pickle(os.path.join(casepath, casedirs["data"], geo_filename), geo_dict)

    geometries_wlegend = {"psPoly":psPoly,"ssPoly":ssPoly,"hk":sortedPoints[ind_hk],"vk":sortedPoints[ind_vk]}
    geometries_nlegend = {"per_y_lower":per_y_lower,"per_y_upper":per_y_upper,}
    plot_geometry_tofile(os.path.join(casepath, casedirs["data"]), geometries_wlegend, geometries_nlegend,
                         "geometry.jpg",zoom=1)

    if verbose:
        plotter = pv.Plotter()
        psPoly = pv.PolyData(np.stack((x_ss, y_ss, np.zeros(len(x_ss)))).T)
        ssPoly = pv.PolyData(np.stack((x_ps, y_ps, np.zeros(len(x_ps)))).T)
        mpslPolyLow = pv.PolyData(np.stack((x_mpsl, y_lower, np.zeros(len(x_mpsl)))).T)
        mpslPolyUpper = pv.PolyData(np.stack((x_mpsl, y_upper, np.zeros(len(x_mpsl)))).T)
        midsPoly = pv.PolyData(np.stack((x_mids, y_mids, np.zeros(len(x_mids)))).T)

        psPoly_asline = lines_from_points(np.stack((x_ss, y_ss, np.zeros(len(x_ss)))).T)
        psPoly_asline["scalars"] = np.arange(psPoly_asline.number_of_points)
        ssPoly_asline = lines_from_points(np.stack((x_ps, y_ps, np.zeros(len(x_ps)))).T)
        ssPoly_asline["scalars"] = np.arange(ssPoly_asline.number_of_points)
        mpslPolyLow_asline = lines_from_points(np.stack((x_mpsl, y_lower, np.zeros(len(x_mpsl)))).T)
        mpslPolyUpper_asline = lines_from_points(np.stack((x_mpsl, y_upper, np.zeros(len(x_mpsl)))).T)
        midsPoly_asline = lines_from_points(np.stack((x_mids, y_mids, np.zeros(len(x_mids)))).T)

        plotter.add_mesh(psPoly, color="red", label="pressure side")
        plotter.add_mesh(ssPoly, color="blue", label="suction side")
        ## plotter.add_mesh(sortedPoints[ind_hk],color="green",point_size=20,label="TE")
        ##plotter.add_mesh(sortedPoints[ind_vk],color="yellow",point_size=20,label="LE")
        plotter.add_mesh(mpslPolyLow)
        plotter.add_mesh(mpslPolyUpper)
        plotter.add_mesh(midsPoly, point_size=5)

        plotter.add_mesh(psPoly_asline)
        plotter.add_mesh(ssPoly_asline)
        plotter.add_mesh(mpslPolyLow_asline)
        plotter.add_mesh(mpslPolyUpper_asline)
        plotter.add_mesh(midsPoly_asline)

        plotter.add_mesh(psPoly, label="psPoly", color="red", point_size=5)
        plotter.add_mesh(ssPoly, label="ssPoly", color="green", point_size=5)
        plotter.add_mesh(psPoly_slice)
        plotter.add_mesh(ssPoly_slice)
        """
        #plotter.add_mesh(pv.PolyData([(-0.004010, 0.048624, 0.000000)]), point_size=20, label="p1")
        #plotter.add_mesh(pv.PolyData([(0.110653, 0.053960, 0.000000)]), point_size=20, label="p2")

        #plotter.add_mesh(pv.PolyData([(0.046967, 0.061204, 0.000000)]), point_size=20, label="p3", color="red")
        plotter.add_mesh(pv.PolyData([(-0.004010, 0.010374, 0.000000)]), point_size=20, label="pt1", color="blue")
        plotter.add_mesh(pv.PolyData(midx_ysouth), point_size=20, label="midx_yss=pt1", color="blue")
        #plotter.add_mesh(pv.PolyData([(0.000733, 0.001322, 0.000000)]), point_size=20, label="pt2")
        #plotter.add_mesh(pv.PolyData([(0.001240, 0.000844, 0.000000)]), point_size=20, label="pt3")
        #plotter.add_mesh(pv.PolyData([(-0.004010, -0.012576, 0.000000)]), point_size=20, label="pt4")
        #plotter.add_mesh(pv.PolyData([(0.110653, 0.069260, 0.000000)]), point_size=20, label="pt5")
        #plotter.add_mesh(pv.PolyData([(0.110653, 0.092210, 0.000000)]), point_size=20, label="pt6")
        #plotter.add_mesh(pv.PolyData([(0.104582, 0.075248, 0.000000)]), point_size=20, label="pt7")
        #plotter.add_mesh(pv.PolyData([(0.104813, 0.074653, 0.000000)]), point_size=20, label="pt8")
        plotter.add_mesh(pv.PolyData([(-0.130000, -0.172116, 0.000000)]), point_size=20, label="pt9")
        plotter.add_mesh(pv.PolyData(inlet_ylow), point_size=20, label="inlet_ylow = pt9")
        #plotter.add_mesh(pv.PolyData([(-0.004010, -0.027876, 0.000000)]), point_size=20, label="pt10")
        plotter.add_mesh(pv.PolyData([(-0.130000, -0.156816, 0.000000)]), point_size=20, label="pt11")
        plotter.add_mesh(pv.PolyData(inlet_southprofile), point_size=20, label="inlet_southprofile = pt11")
        plotter.add_mesh(pv.PolyData([(-0.130000, -0.133866, 0.000000)]), point_size=20, label="pt12")
        plotter.add_mesh(pv.PolyData(inlet_northprofile), point_size=20, label="inlet_northprofile = pt12")
        plotter.add_mesh(pv.PolyData([(-0.130000, -0.095616, 0.000000)]), point_size=20, label="pt13")
        plotter.add_mesh(pv.PolyData(inlet_yhigh), point_size=20, label="inlet_high = pt13")
        #plotter.add_mesh(pv.PolyData([(-0.004010, -0.027876, 0.000000)]), point_size=20, label="pt14")
        #plotter.add_mesh(pv.PolyData([(0.110653, 0.053960, 0.000000)]), point_size=20, label="pt15")
        #plotter.add_mesh(pv.PolyData([(0.110653, 0.053960, 0.000000)]), point_size=20, label="pt16")
        #plotter.add_mesh(pv.PolyData([(0.260000, 0.110953, 0.000000)]), point_size=20, label="pt17")
        plotter.add_mesh(pv.PolyData([(0.260000, 0.126253, 0.000000)]), point_size=20, label="pt18")
        plotter.add_mesh(pv.PolyData(outlet_southprofile), point_size=20, label="outlet_southprofile=pt18")
        plotter.add_mesh(pv.PolyData([(0.260000, 0.149203, 0.000000)]), point_size=20, label="pt19")
        plotter.add_mesh(pv.PolyData(outlet_northprofile), point_size=20, label="outlet_northprofile=pt19")
        #plotter.add_mesh(pv.PolyData([(0.260000, 0.187453, 0.000000)]), point_size=20, label="pt20")
        plotter.add_mesh(pv.PolyData([(0.110653, 0.130460, 0.000000)]), point_size=20, label="pt21")
        plotter.add_mesh(pv.PolyData(outlet_yhigh), point_size=20, label="outlet_high = pt21")
        """
        plotter.show_axes()
        plotter.add_legend()
        plotter.show()

    return geo_dict


def create_geometry_fromnacaairfoil(nacadigits, numberofpoints, finite_TE, half_cosine_spacing, x_inlet, x_outlet,
                                    pitch, geoscaling, blade_shift, staggerangle, span_z, casepath, verbose=False):
    # =============================================================================
    # Bestimmung Profilparameter
    # =============================================================================
    sortedPoints, psPoly, ssPoly, ind_vk, ind_hk, midsPoly, beta_meta_01, beta_meta_02, camberangle = create_naca_geoparas(
        nacadigits,
        numberofpoints,
        finite_TE,
        half_cosine_spacing,
        geoscaling,
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

    writeTecplot1DFile(os.path.join(casepath, '01_Meshing', 'geom.dat'), ['x', 'z'],
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
    write_pickle(os.path.join(casepath, casedirs["data"], geo_filename), geo_dict)

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


def create_naca_geoparas(nacadigits, numberofpoints, finite_TE, half_cosine_spacing, geoscaling, staggerangle,
                         verbose=True):
    ptsx, ptsy = naca(nacadigits, numberofpoints, finite_TE, half_cosine_spacing)
    ind_hk = 0
    ind_vk = numberofpoints
    points = np.stack((ptsx[:-1], ptsy[:-1], np.zeros(numberofpoints * 2) - 1)).T
    scalegeo = geoscaling / vecAbs(points[ind_hk] - points[ind_vk])
    poly = pv.PolyData(points)
    poly.points *= scalegeo
    poly.rotate_z(staggerangle)
    points = poly.points
    ssPoly, psPoly = extractSidePolys(ind_hk, ind_vk, poly, verbose)
    midsPoly = midline_from_sides(ind_hk, ind_vk, points, psPoly, ssPoly)
    metal_angle_hk, metal_angle_vk, camber_angle = angles_from_mids(midsPoly)

    if verbose:
        p = pv.Plotter()
        p.add_mesh(points, color="orange", label="points")
        p.add_mesh(psPoly, color="green", label="psPoly")
        p.add_mesh(ssPoly, color="black", label="ssPoly")
        p.add_mesh(midsPoly, color="black", label="midsPoly")
        p.add_legend()
        p.show()
    return points, psPoly, ssPoly, ind_vk, ind_hk, midsPoly, metal_angle_vk, metal_angle_hk, camber_angle


def extract_geo_paras(points, alpha, verbose=False):
    """
    This function is extracting profile-data as stagger-angle, midline, psPoly, ssPoly and more from a set of points
    Be careful, you need a suitable alpha-parameter in order to get the right geometry
    The calculation of the leading-edge and trailing-edge index needs time and its not 100% reliable (yet)
    Keep in mind, to check the results!
    :param points: array of points in 3d with the shape (n,3)
    :param alpha: nondimensional alpha-coefficient (calcConcaveHull)
    :param verbose: bool for plots
    :return: points, psPoly, ssPoly, ind_vk, ind_hk, midsPoly, metal_angle_vk, metal_angle_hk
    """

    origPoly = pv.PolyData(points)
    xs, ys = calcConcaveHull(points[:, 0], points[:, 1], alpha)
    points = np.stack((xs, ys, np.zeros(len(xs)))).T
    sortedPoly = pv.PolyData(points)

    ind_hk, ind_vk, veronoi_mid = extract_vk_hk(origPoly, sortedPoly)
    psPoly, ssPoly = extractSidePolys(ind_hk, ind_vk, sortedPoly, verbose)
    midsPoly = midline_from_sides(ind_hk, ind_vk, points, psPoly, ssPoly)
    metal_angle_hk, metal_angle_vk, camber_angle = angles_from_mids(midsPoly)

    if verbose:
        p = pv.Plotter()
        p.add_mesh(points, color="orange", label="points")
        p.add_mesh(psPoly, color="green", label="psPoly")
        p.add_mesh(ssPoly, color="black", label="ssPoly")
        p.add_mesh(midsPoly, color="black", label="midsPoly")
        p.add_mesh(veronoi_mid, color="yellow", label="veronoi_mid")
        p.add_legend()
        p.show()

    return points, psPoly, ssPoly, ind_vk, ind_hk, midsPoly, metal_angle_vk, metal_angle_hk, camber_angle


def run_create_geometry(settings_yaml):
    case_path = os.path.abspath(os.path.dirname(settings_yaml))
    settings = yaml_dict_read(settings_yaml)

    meshpath = os.path.join(case_path, casedirs["meshing"])
    datpath = os.path.join(case_path, casedirs["data"])
    if not os.path.isdir(meshpath):
        os.mkdir(meshpath)
    if not os.path.isdir(datpath):
        os.mkdir(datpath)

    print("create_geometry")
    if settings["geometry"]["algorithm"] == "from_pointcloud":
        ptstxtfile = os.path.join(os.path.abspath(case_path), settings["geometry"]["ptcloud_profile"])

        geo_dict = create_geometry_frompointcloud(ptstxtfile,
                                                  settings, case_path)

    elif settings["geometry"]["algorithm"] == "naca_airfoil_generator":
        geo_dict = create_geometry_fromnacaairfoil(settings["geometry"]["naca_digits"],
                                                   settings["geometry"]["numberofpoints"],
                                                   settings["geometry"]["finite_TE"],
                                                   settings["geometry"]["half_cosine_spacing"],
                                                   settings["geometry"]["x_inlet"],
                                                   settings["geometry"]["x_outlet"],
                                                   settings["geometry"]["pitch"],
                                                   settings["geometry"]["camberlength"],
                                                   settings["geometry"]["shift_domain"],
                                                   settings["geometry"]["staggerangle"],
                                                   settings["mesh"]["extrudeLength"],
                                                   case_path, )

    elif settings["geometry"]["algorithm"] == "prep":
        geo_dict = prep_geo(settings_yaml)

    # blockpoints(geo_dict,settings)
    return 0
