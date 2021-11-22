import numpy as np
import pyvista as pv
from matplotlib import path as mpltPath
from scipy.interpolate import splprep, splev
from scipy.spatial import Voronoi
from itertools import product

from NTR.utils.geom_functions.pointcloud import calcConcaveHull
from NTR.utils.geom_functions.spline import refine_spline, calcMidPoints
from NTR.utils.geom_functions.distance import closest_node_index, distant_node_index
from NTR.utils.mathfunctions import vecAbs, angle_between, vecDir
from NTR.utils.mesh_handling.pyvista_utils import lines_from_points, polyline_from_points
from NTR.utils.geom_functions.distance import calc_largedistant_idx


def extract_vk_hk(origPoly, sortedPoly, verbose=False):
    """
    This function is calculating the leading-edge and trailing edge of a long 2d-body
    The function is not 100% reliable yet. The computation is iterative and it can take a while
    Points in origPoly and sortedPoly have to have defined points on the LE and TE, otherwise a LE or TE is not defined
    and it will be random which point will be found near the LE / TE
    :param origPoly: all original points, unsorted
    :param sortedPoly: sorted via calcConcaveHull
    :param verbose: bool (True -> plots, False -> silent)
    :return: returns indexes of LE(vk) and TE(hk) from sortedPoints
    """

    xs, ys = sortedPoly.points[::, 0], sortedPoly.points[::, 1]
    ind_1, ind_2 = calc_largedistant_idx(xs, ys)

    nopt = sortedPoly.number_of_points

    allowed_shift = int(nopt * 0.08)
    shifts = np.arange(-allowed_shift,allowed_shift+1)
    ind_1_ts = (shifts + ind_1)%nopt
    ind_2_ts = (shifts + ind_2)%nopt

    combs = list(product(ind_1_ts,ind_2_ts))

    def midLength(ind_1, ind_2 ,verbose=False):
        psPoly, ssPoly = extractSidePolys(ind_1, ind_2, sortedPoly, False)
        midsPoly = midline_from_sides(ind_1, ind_2, sortedPoly.points, psPoly, ssPoly)
        arclength = midsPoly.compute_arc_length()["arc_length"]
        midslength = sum(arclength)

        if verbose:
            p = pv.Plotter()
            p.add_mesh(midsPoly)
            p.add_mesh(sortedPoly)
            p.add_mesh(sortedPoly.points[ind_1], color="yellow", point_size=20)
            p.add_mesh(sortedPoly.points[ind_2], color="red", point_size=20)
            p.add_text("length: " + str(midslength))

            p.show()
        return midslength

    midLengths = []
    for ind_1_t, ind2_t in combs:
        midLengths.append(midLength(ind_1_t, ind2_t,verbose))
    cids = midLengths.index(max(midLengths))
    ind_1, ind_2 = combs[cids]

    if sortedPoly.points[ind_1][0] > sortedPoly.points[ind_2][0]:
        ind_vk = ind_2
        ind_hk = ind_1
    else:
        ind_vk = ind_1
        ind_hk = ind_2
    return ind_hk, ind_vk


def extract_edge_poi(try_center, try_radius, mids, direction, sortedPoly, verbose=False):
    mids_minx = mids[[i[0] for i in mids].index(min([i[0] for i in mids]))]
    mids_maxx = mids[[i[0] for i in mids].index(max([i[0] for i in mids]))]

    mids_tangent = mids_minx - mids_maxx

    splitBoxLength = vecAbs(try_center - sortedPoly.points[distant_node_index(try_center, sortedPoly.points)]) * 2.1
    splitBox = pv.Plane(center=(0, 0, 0), direction=(0, 0, 1), i_size=try_radius * 1.6, j_size=splitBoxLength,
                        i_resolution=100, j_resolution=100)

    rotate = -angle_between(mids_tangent, np.array([0, 1, 0])) / np.pi * 180
    if direction == "low":
        splitBox = pv.PolyData(np.array([i for i in splitBox.points if i[1] <= 0]))
    elif direction == "high":
        splitBox = pv.PolyData(np.array([i for i in splitBox.points if i[1] >= 0]))

    splitBox = splitBox.delaunay_2d()
    splitBox = splitBox.extrude((0, 0, 0.1))
    splitBox.translate((0, 0, -0.05))

    if direction == "low":
        splitBox.rotate_z(-rotate)
        splitBox.translate(vecDir(try_center - mids_minx) * try_radius)
    elif direction == "high":
        splitBox.rotate_z(-rotate)
        splitBox.translate(vecDir(try_center - mids_maxx) * try_radius)

    splitBox.points += try_center
    enclosedBoxPoints = sortedPoly.select_enclosed_points(splitBox)
    checkPoints = [i for idx, i in enumerate(enclosedBoxPoints.points) if
                   enclosedBoxPoints["SelectedPoints"][idx] == 1]

    if verbose:
        p = pv.Plotter()
        p.add_mesh(sortedPoly)
        p.add_mesh(pv.PolyData(np.asarray(checkPoints)), color="blue")
        p.add_mesh(splitBox.extract_feature_edges())
        p.add_mesh(np.array(mids), color="red")
        p.show()

    return checkPoints


def extractSidePolys(ind_hk, ind_vk, sortedPoly, verbose=False):
    xs, ys = list(sortedPoly.points[::, 0]), list(sortedPoly.points[::, 1])

    if ind_vk < ind_hk:
        x_ss = xs[ind_vk:ind_hk + 1]
        y_ss = ys[ind_vk:ind_hk + 1]

        y_ps = ys[ind_hk:] + ys[:ind_vk + 1]
        x_ps = xs[ind_hk:] + xs[:ind_vk + 1]

    else:
        x_ss = xs[ind_hk:ind_vk + 1]
        y_ss = ys[ind_hk:ind_vk + 1]

        y_ps = ys[ind_vk:] + ys[:ind_hk + 1]
        x_ps = xs[ind_vk:] + xs[:ind_hk + 1]

    psl_helper = polyline_from_points(np.stack((x_ps, y_ps, np.zeros(len(x_ps)))).T)
    ssl_helper = polyline_from_points(np.stack((x_ss, y_ss, np.zeros(len(x_ss)))).T)

    if psl_helper.length > ssl_helper.length:

        psPoly = pv.PolyData(ssl_helper.points)
        ssPoly = pv.PolyData(psl_helper.points)
    else:

        psPoly = pv.PolyData(psl_helper.points)
        ssPoly = pv.PolyData(ssl_helper.points)
    if verbose:
        p = pv.Plotter()
        psl = polyline_from_points(psPoly.points)
        psl.point_arrays["scalars"] = range(psl.number_of_points)
        p.add_mesh(psl, label="psPoly")
        ssl = polyline_from_points(ssPoly.points)
        ssl.point_arrays["scalars"] = range(ssl.number_of_points)
        p.add_mesh(ssl, label="ssPoly")
        p.add_mesh(psPoly, label="psPoly", color="white", point_size=1)
        p.add_mesh(ssPoly, label="ssPoly", color="black", point_size=1)
        p.add_legend()
        p.add_axes()
        p.show()

    return ssPoly, psPoly


def calcMidPassageStreamLine(x_mcl, y_mcl, beta1, beta2, x_inlet, x_outlet, t):
    """
    Returns mid-passage line from sceletal-line
    Returns two lists of Points representing a curve through the passage


    Input:
    x_mcl, y_mcl = Tuple
    beta1, beta2 = Angle in deg - Beta = Anströmwinkel
    x_inlet, x_outlet = scalar - representing position x-component of in/outlet
    t = scalar pitch
    """

    delta_x_vk = x_mcl[0] - x_inlet
    delta_y_vk = np.tan(np.deg2rad(beta1 - 90)) * delta_x_vk

    p_inlet_x = x_mcl[0] - delta_x_vk
    p_inlet_y = y_mcl[0] - delta_y_vk

    delta_x_hk = x_outlet - x_mcl[-1]
    delta_y_hk = delta_x_hk * np.tan(np.deg2rad(beta2 - 90))

    p_outlet_x = x_mcl[-1] + delta_x_hk
    p_outlet_y = y_mcl[-1] + delta_y_hk

    x_mpsl = [p_inlet_x] + list(x_mcl) + [p_outlet_x]
    y_mpsl = [p_inlet_y] + list(y_mcl) + [p_outlet_y]

    for i in range(len(x_mpsl)):
        y_mpsl[i] = y_mpsl[i] + 0.5 * t

    return refine_spline(x_mpsl, y_mpsl, 1000)


def sortProfilePoints(x, y, alpha):
    x, y = calcConcaveHull(x, y, alpha=alpha)

    ind_vk = x.index(min(x))

    x_by_vk = None
    y_by_vk = None

    if ind_vk != 0:

        if y[ind_vk] < y[ind_vk - 1]:
            x_by_vk = x[ind_vk:] + x[:ind_vk + 1]
            y_by_vk = y[ind_vk:] + y[:ind_vk + 1]

        else:
            x_by_vk = x[ind_vk:] + x[:ind_vk + 1]
            y_by_vk = y[ind_vk:] + y[:ind_vk + 1]

    else:

        if y[ind_vk] > y[ind_vk + 1]:
            x_by_vk = x_by_vk[::-1]
            y_by_vk = y_by_vk[::-1]

    ind_hk = x_by_vk.index(max(x_by_vk))

    x_ss = []
    y_ss = []

    x_ps = []
    y_ps = []

    for i in range(len(x_by_vk)):

        if i <= ind_hk:
            x_ss.append(x_by_vk[i])
            y_ss.append(y_by_vk[i])
        else:
            x_ps.append(x_by_vk[i])
            y_ps.append(y_by_vk[i])

    if y_ss[1] < y_ss[0]:
        x = list(x_ss)
        y = list(y_ss)

        x_ss = list(x_ps)
        y_ss = list(y_ps)

        x_ps = x
        y_ps = y

    return x_ss, y_ss, x_ps, y_ps


def calc_vk_hk(x_koords, y_koords, beta_01, beta_02):
    # Vorderkante

    # Punkt bestimmen deutlich unterhalt der Vorderkante liegt
    p01 = [min(x_koords) - (max(x_koords) - min(x_koords)), min(y_koords) - (max(x_koords) - min(x_koords))]

    # Vektor bestimmen der Richtung des vorgegebenen Winkels der Vorderkante besitzt
    vec_01 = [np.sin(np.deg2rad(180.0 - beta_01)) * 0.01, np.cos(np.deg2rad(180.0 - beta_01)) * 0.01]

    # Vektor bestimmen der orthogonal dazu ist
    vec_02 = [1, -vec_01[0] * 1 / vec_01[1]]

    # alle Koords durchgehen
    dists = []

    betrag_vec_02 = np.sqrt(vec_02[0] ** 2 + vec_02[1] ** 2)
    for i in range(len(y_koords)):
        delta_x = x_koords[i] - p01[0]
        delta_y = y_koords[i] - p01[1]
        d = np.sqrt((vec_02[0] * delta_y - delta_x * vec_02[1]) ** 2) / betrag_vec_02
        dists.append(d)

    index_vk = np.argmin(dists)

    # Hinterkante
    p01 = [max(x_koords) + (max(x_koords) - min(x_koords)), min(y_koords) - (max(x_koords) - min(x_koords))]
    vec_01 = [-np.sin(np.deg2rad(beta_02)) * 0.01, np.cos(np.deg2rad(beta_02)) * 0.01]
    vec_02 = [1, -vec_01[0] * 1 / vec_01[1]]
    dists = []

    betrag_vec_02 = np.sqrt(vec_02[0] ** 2 + vec_02[1] ** 2)
    for i in range(len(y_koords)):
        delta_x = x_koords[i] - p01[0]
        delta_y = y_koords[i] - p01[1]
        d = np.sqrt((vec_02[0] * delta_y - delta_x * vec_02[1]) ** 2) / betrag_vec_02
        dists.append(d)

    index_hk = np.argmin(dists)

    return index_vk, index_hk


def midline_from_sides(ind_hk, ind_vk, points, psPoly, ssPoly):
    x_ps, y_ps = psPoly.points[::, 0], psPoly.points[::, 1]
    x_ss, y_ss = ssPoly.points[::, 0], ssPoly.points[::, 1]

    midsres = 100
    if x_ps[0] < x_ps[-1]:
        ax, ay = refine_spline(x_ps[::-1], y_ps[::-1], midsres)
    else:
        ax, ay = refine_spline(x_ps, y_ps, midsres)
    if x_ss[0] < x_ss[-1]:
        bx, by = refine_spline(x_ss[::-1], y_ss[::-1], midsres)
    else:
        bx, by = refine_spline(x_ss, y_ss, midsres)
    xmids, ymids = ((ax + bx) / 2, (ay + by) / 2)
    xmids = np.array(xmids)[::-1][1:-1]
    ymids = np.array(ymids)[::-1][1:-1]
    xmids[0] = points[ind_vk][0]
    ymids[0] = points[ind_vk][1]
    xmids[-1] = points[ind_hk][0]
    ymids[-1] = points[ind_hk][1]
    midsPoly = lines_from_points(np.stack((xmids, ymids, np.zeros(len(ymids)))).T)
    return midsPoly


def angles_from_mids(midsPoly):
    xmids, ymids = midsPoly.points[::, 0], midsPoly.points[::, 1]
    vk_tangent = np.stack((xmids[0] - xmids[1], ymids[0] - ymids[1], 0)).T
    hk_tangent = np.stack((xmids[-2] - xmids[-1], ymids[-2] - ymids[-1], 0)).T
    camber = np.stack((xmids[0] - xmids[-1], ymids[0] - ymids[-1], 0)).T
    metal_angle_vk = angle_between(vk_tangent, np.array([0, 1, 0])) / np.pi * 180
    metal_angle_hk = angle_between(hk_tangent, np.array([0, 1, 0])) / np.pi * 180
    camber_angle = angle_between(camber, np.array([0, 1, 0])) / np.pi * 180
    return metal_angle_hk, metal_angle_vk, camber_angle


def calcMeanCamberLine(x, y, alpha):
    # vk und hk bestimmen

    ind_vk, ind_hk = calc_vk_hk(x, y)

    x_vk = x[ind_vk]
    y_vk = y[ind_vk]

    x_hk = x[ind_hk]
    y_hk = y[ind_hk]

    x_ss, y_ss, x_ps, y_ps = sortProfilePoints(x, y, alpha=alpha)

    x_mid_ss, y_mid_ss = calcMidPoints(x_ss, y_ss, x_ps, y_ps)
    x_mid_ps, y_mid_ps = calcMidPoints(x_ps, y_ps, x_ss, y_ss)
    x_mids, y_mids = calcMidPoints(x_mid_ps, y_mid_ps, x_mid_ss, y_mid_ss)

    return x_mids, y_mids, x_ss, y_ss, x_ps, y_ps, x_vk, y_vk, x_hk, y_hk


def extract_profile_paras(points, verbose=False):
    points2d = points[:, 0:2]
    vor = Voronoi(points2d)

    midline = []
    path = mpltPath.Path(points2d)
    for idx, r in enumerate(vor.regions):

        pts = vor.vertices[r]
        pts3d = np.insert(pts, 2, 0, axis=1)
        inside = path.contains_points(pts)
        pts3dclean = [pts3d[idx] for idx, i in enumerate(pts3d) if inside[idx] == True]
        for p in pts3dclean:
            if not p[0] in [i[0] for i in midline]:
                midline.append(p)

    midpoints = pv.PolyData(midline)

    xsortedpoints = midpoints.points[np.argsort(midpoints.points[:, 0])]

    twodpts = xsortedpoints[:, 0:2].T

    (tck, u), fp, ier, msg = splprep(twodpts, u=None, per=0, k=5, full_output=True)
    # s = optional parameter (default used here)
    # #print('Spline score:',fp) #goodness of fit flatlines after a given s value (and higher), which captures the default s-value as well
    x_new, y_new = splev(u, tck, der=0)

    splineNew = np.stack((x_new, y_new)).T  # np.zeros(len(x_new)

    inside = path.contains_points(splineNew)
    splinePointsNewClean = [splineNew[idx] for idx, i in enumerate(splineNew) if inside[idx] == True]
    splineNewClean = np.insert(splinePointsNewClean, 2, 0, axis=1)

    firstBase = splineNewClean[0]
    firstBaseUp = splineNewClean[1]

    secondBase = splineNewClean[-1]
    secondBaseUp = splineNewClean[-2]

    firstDir = vecDir(firstBase - firstBaseUp)
    secondDir = vecDir(secondBase - secondBaseUp)

    extend_first = pv.Line(firstBase, firstBase + firstDir, 100000)
    extend_second = pv.Line(secondBase, secondBaseUp + secondDir, 100000)

    bladeLine = polyline_from_points(np.vstack([points, [points[0]]]))

    firsthitpoint = extend_first.slice_along_line(bladeLine).points[0]
    secondhitpoint = extend_second.slice_along_line(bladeLine).points[0]

    hitids = [closest_node_index(firsthitpoint, bladeLine.points), closest_node_index(secondhitpoint, bladeLine.points)]

    low_hitid = min(hitids)
    high_hitid = max(hitids)

    hitpoint_lowid = secondhitpoint  # bladeLine.points[low_hitid]
    hitpoint_highid = firsthitpoint  # bladeLine.points[high_hitid]

    if bladeLine.points[low_hitid][0] < hitpoint_highid[0]:
        vk_point = hitpoint_lowid  # bladeLine.points[vk_point_id]
        hk_dir = extend_first

        hk_point = hitpoint_highid  # bladeLine.points[hk_point_id]
        vk_dir = extend_second
    else:
        hk_point = hitpoint_lowid  # bladeLine.points[vk_point_id]
        vk_dir = extend_first
        vk_point = hitpoint_highid  # bladeLine.points[hk_point_id]
        hk_dir = extend_second

    ss_poly = pv.PolyData(bladeLine.points[low_hitid:high_hitid])
    ps_poly = pv.PolyData(np.append(bladeLine.points[:low_hitid], bladeLine.points[high_hitid:], axis=0))

    beta_01 = angle_between(vk_dir.points[0] - vk_dir.points[-1], np.array([1, 0, 0])) * 360 / (2 * np.pi)
    beta_02 = angle_between(hk_dir.points[0] - hk_dir.points[-1], np.array([1, 0, 0])) * 360 / (2 * np.pi)

    centerline = lines_from_points(
        np.append(np.append(np.asarray([vk_point]), splineNewClean, axis=0), np.asarray([hk_point]), axis=0))

    if verbose:
        plotter = pv.Plotter()
        plotter.add_mesh(ss_poly, color="blue")
        plotter.add_mesh(ps_poly, color="yellow")
        plotter.add_mesh(pv.PolyData(vk_point), point_size=20, color="orange")
        plotter.add_mesh(pv.PolyData(hk_point), point_size=20, color="grey")

        plotter.add_mesh(centerline, point_size=20)
        plotter.show()

    return ss_poly.points, ps_poly.points, centerline.points, beta_01, beta_02
