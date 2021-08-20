import math
import numpy as np
import pyvista as pv
from matplotlib import path as mpltPath

from scipy.interpolate import splprep, interp1d, splev
from scipy.spatial import Delaunay, distance
from scipy.optimize import minimize
from scipy.spatial.qhull import Voronoi
from scipy.spatial.distance import squareform, pdist

from NTR.utils.mathfunctions import vecAbs, vecDir, closest_node_index, angle_between, splineCurvature
from NTR.utils.thermoFunctions import Sutherland_Law
from NTR.utils.boundaryLayerFunctions import calcWallShearStress
from NTR.utils.simFunctions import sort_values_by_pitch
from NTR.utils.pyvista_utils import slice_midspan_z, polyline_from_points, lines_from_points


def refine_spline(x, y, res):
    """
    https://stackoverflow.com/questions/51512197/python-equidistant-points-along-a-line-joining-set-of-points/51515357

    :param x:
    :param y:
    :param res:
    :return: x,y
    """

    distance = np.cumsum(np.sqrt(np.ediff1d(x, to_begin=0) ** 2 + np.ediff1d(y, to_begin=0) ** 2))
    distance = distance / distance[-1]

    fx, fy = interp1d(distance, x), interp1d(distance, y)

    alpha_ = np.linspace(0, 1, res)
    x, y = fx(alpha_), fy(alpha_)
    return x, y


def midpoint(x1, y1, x2, y2):
    return ((x1 + x2) / 2, (y1 + y2) / 2)


def calcMidPoints(x1, y1, x2, y2, res=100, helperpoints=50000):
    # The reason is the tolerance which needs to be low enough for a reasonable realiability
    # but on the same side this is an approximation and the tolerance cant be too low
    res = res
    res_oppposit = helperpoints

    x1_one, y1_one = refine_spline(x1, y1, res)
    x2_one, y2_one = refine_spline(x2, y2, res_oppposit)

    x_mid_ss = []
    y_mid_ss = []

    for i in range(len(x1_one)):

        dists = []

        for j in range(len(x2_one)):
            dist = ((x1_one[i] - x2_one[j]) ** 2 + (y1_one[i] - y2_one[j]) ** 2) ** (0.5)

            dists.append(dist)

        index_p = np.argmin(dists)

        p_x = x2_one[index_p]
        p_y = y2_one[index_p]

        x_mid, y_mid = midpoint(p_x, p_y, x1_one[i], y1_one[i])

        dists_ss = []
        dists_ps = []

        for j in range(len(x2_one)):
            dists_ss.append(((x_mid - x2_one[j]) ** 2 + (y_mid - y2_one[j]) ** 2) ** (0.5))
        for j in range(len(x1_one)):
            dists_ps.append(((x_mid - x1_one[j]) ** 2 + (y_mid - y1_one[j]) ** 2) ** (0.5))

        dist_1 = dists_ss[np.argmin(dists_ss)]
        dist_2 = dists_ps[np.argmin(dists_ps)]
        if abs(dist_1-dist_2)/(dist_1+dist_2) < 0.05:#, rtol=tolerance
            x_mid_ss.append(x_mid)
            y_mid_ss.append(y_mid)

    x_mid_ss, y_mid_ss = zip(*sorted(zip(x_mid_ss, y_mid_ss)))

    return np.array(x_mid_ss), np.array(y_mid_ss)


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


def calcConcaveHull_optimize(xs, ys, alpha):
    def func(alpha):
        xss, yss = calcConcaveHull(xs, ys, alpha)
        return len(xs) - len(xss)

    res = minimize(func, alpha, method='Powell',
                   options={'xatol': 1e-8, 'disp': True})

    opt_x, opt_y = calcConcaveHull(xs, ys, res.x[0])
    return opt_x, opt_y


def calcConcaveHull(x, y, alpha):
    """
    origin: https://stackoverflow.com/questions/50549128/boundary-enclosing-a-given-set-of-points/50714300#50714300
    """
    points = []
    for i in range(len(x)):
        points.append([x[i], y[i]])

    points = np.asarray(points)

    def alpha_shape(points, alpha, only_outer=True):
        """
        Compute the alpha shape (concave hull) of a set of points.
        :param points: np.array of shape (n,2) points.
        :param alpha: alpha value.
        :param only_outer: boolean value to specify if we keep only the outer border
        or also inner edges.
        :return: set of (i,j) pairs representing edges of the alpha-shape. (i,j) are
        the indices in the points array.
        """

        assert points.shape[0] > 3, "Need at least four points"

        def add_edge(edges, i, j):
            """
            Add an edge between the i-th and j-th points,
            if not in the list already
            """
            if (i, j) in edges or (j, i) in edges:
                # already added
                assert (j, i) in edges, "Can't go twice over same directed edge right?"
                if only_outer:
                    # if both neighboring triangles are in shape, it's not a boundary edge
                    edges.remove((j, i))
                return
            edges.add((i, j))

        tri = Delaunay(points)
        edges = set()
        # Loop over triangles:
        # ia, ib, ic = indices of corner points of the triangle
        for ia, ib, ic in tri.vertices:
            pa = points[ia]
            pb = points[ib]
            pc = points[ic]
            # Computing radius of triangle circumcircle
            # www.mathalino.com/reviewer/derivation-of-formulas/derivation-of-formula-for-radius-of-circumcircle
            a = np.sqrt((pa[0] - pb[0]) ** 2 + (pa[1] - pb[1]) ** 2)
            b = np.sqrt((pb[0] - pc[0]) ** 2 + (pb[1] - pc[1]) ** 2)
            c = np.sqrt((pc[0] - pa[0]) ** 2 + (pc[1] - pa[1]) ** 2)
            s = (a + b + c) / 2.0

            area = np.sqrt(s * (s - a) * (s - b) * (s - c))

            circum_r = a * b * c / (4.0 * area)
            if circum_r < alpha:
                add_edge(edges, ia, ib)
                add_edge(edges, ib, ic)
                add_edge(edges, ic, ia)
        return edges

    def find_edges_with(i, edge_set):
        i_first = [j for (x, j) in edge_set if x == i]
        i_second = [j for (j, x) in edge_set if x == i]
        return i_first, i_second

    def stitch_boundaries(edges):
        edge_set = edges.copy()
        boundary_lst = []
        while len(edge_set) > 0:
            boundary = []
            edge0 = edge_set.pop()
            boundary.append(edge0)
            last_edge = edge0
            while len(edge_set) > 0:
                i, j = last_edge
                j_first, j_second = find_edges_with(j, edge_set)
                if j_first:
                    edge_set.remove((j, j_first[0]))
                    edge_with_j = (j, j_first[0])
                    boundary.append(edge_with_j)
                    last_edge = edge_with_j
                elif j_second:
                    edge_set.remove((j_second[0], j))
                    edge_with_j = (j, j_second[0])  # flip edge rep
                    boundary.append(edge_with_j)
                    last_edge = edge_with_j

                if edge0[0] == last_edge[1]:
                    break

            boundary_lst.append(boundary)
        return boundary_lst

    edges = alpha_shape(points, alpha)
    boundary_lst = stitch_boundaries(edges)
    x_new = []
    y_new = []

    for i in range(len(boundary_lst[0])):
        x_new.append(points[boundary_lst[0][i][0]][0])
        y_new.append(points[boundary_lst[0][i][0]][1])

    return x_new, y_new


def veronoi_midline(points, verbose=True):
    points2d = points[::, 0:2]
    vor = Voronoi(points2d)
    midline = []
    for idx, r in enumerate(vor.regions):

        pts = vor.vertices[r]
        pts3d = np.insert(pts, 2, 0, axis=1)

        inside = inside_poly(points2d, pts3d[::, 0:2])
        pts3dclean = [i for idx, i in enumerate(pts3d) if inside[idx] == True]
        for p in pts3dclean:
            if not p[0] in [i[0] for i in midline]:
                midline.append(p)

    midpoints = pv.PolyData(midline)

    xsortedpoints = midpoints.points[np.argsort(midpoints.points[:, 0])]

    twodpts = xsortedpoints[:, 0:2].T

    (tck, u), fp, ier, msg = splprep(twodpts, u=None, per=0, k=5,s=100, full_output=True)

    x_new, y_new = splev(u, tck, der=0)

    x_new, y_new = refine_spline(x_new, y_new, 100)
    splineNew = np.stack((x_new, y_new, np.zeros(len(x_new)))).T

    inside = inside_poly(points2d, splineNew[::, 0:2])
    splineNewclean = [i for idx, i in enumerate(splineNew) if inside[idx] == True]
    splines = []
    for p in splineNewclean:
        if not p[0] in [i[0] for i in splines]:
            splines.append(p)

    splines = polyline_from_points(np.array(splines))

    while max(splineCurvature(splines.points[:,0],splines.points[:,1])>max(splineCurvature(points2d[:,0],points2d[:,1]))):
        splines.point_arrays["curvature"] = splineCurvature(splines.points[:,0],splines.points[:,1])
        delid = np.where(splines["curvature"]>max(splineCurvature(points2d[:,0],points2d[:,1])))[0]
        pts = np.asarray([p for idp, p in enumerate(splines.points) if idp not in delid])
        splines = polyline_from_points(pts)

    if verbose:
        p = pv.Plotter()
        p.add_mesh(points)
        p.add_mesh(splines)
        p.show()
    return splines


def line_intersection(point_a1, point_a2,
                      point_b1, point_b2):
    def det_2d(a, b):
        return a[0] * b[1] - a[1] * b[0]

    xdiff = (point_a1[0] - point_a2[0], point_b1[0] - point_b2[0])
    ydiff = (point_a1[1] - point_a2[1], point_b1[1] - point_b2[1])

    div = det_2d(xdiff, ydiff)
    if div == 0:
        return None

    d = (det_2d(point_a1, point_a2), det_2d(point_b1, point_b2))
    x = det_2d(d, xdiff) / div
    y = det_2d(d, ydiff) / div
    return x, y


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


def rotate_points(origin, x, y, angle):
    """
    Rotate points counterclockwise (rotation around z-axis) by a given angle around a given origin.
    """

    angle = angle * (math.pi / 180.0)

    ox, oy = origin

    new_x = []
    new_y = []

    for i in range(len(x)):
        qx = ox + math.cos(angle) * (x[i] - ox) - math.sin(angle) * (y[i] - oy)
        qy = oy + math.sin(angle) * (x[i] - ox) + math.cos(angle) * (y[i] - oy)

        new_x.append(qx)
        new_y.append(qy)

    return new_x, new_y


def calc_largedistant_idx(x_koords, y_koords):
    A = np.dstack((x_koords, y_koords))[0]
    D = squareform(pdist(A))
    #    N = np.max(D)
    I = np.argmax(D)
    I_row, I_col = np.unravel_index(I, D.shape)

    index_1 = I_row
    index_2 = I_col

    return index_1, index_2


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
    def extract_edge_poi(try_center, try_radius, mids, direction, sortedPoly, verbose=False):
        mids_minx = mids[[i[0] for i in mids].index(min([i[0] for i in mids]))]
        mids_maxx = mids[[i[0] for i in mids].index(max([i[0] for i in mids]))]

        mids_tangent = mids_minx - mids_maxx

        splitBoxLength = vecAbs(try_center - sortedPoly.points[distant_node_index(try_center, sortedPoly.points)])*2.1
        splitBox = pv.Plane(center=(0, 0, 0), direction=(0, 0, 1), i_size=try_radius*1.6, j_size=splitBoxLength,
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
            splitBox.translate(mids_maxx-mids_minx)
        elif direction == "high":
            splitBox.rotate_z(-rotate)
            splitBox.translate(mids_minx-mids_maxx)

        splitBox.points += try_center
        enclosedBoxPoints = sortedPoly.select_enclosed_points(splitBox)
        checkPoints = [i for idx, i in enumerate(enclosedBoxPoints.points) if
                       enclosedBoxPoints["SelectedPoints"][idx] == 1]

        if verbose:
            p=pv.Plotter()
            p.add_mesh(sortedPoly)
            p.add_mesh(pv.PolyData(np.asarray(checkPoints)),color="blue")
            p.add_mesh(splitBox.extract_feature_edges())
            p.add_mesh(np.array(mids),color="red")
            p.show()
        return checkPoints

    xs, ys = sortedPoly.points[::,0], sortedPoly.points[::,1]
    x_new, y_new = refine_spline(xs, ys, 10000)
    splineNew = np.stack((x_new, y_new, np.zeros(len(x_new)))).T
    linePoly = lines_from_points(splineNew)
    veronoi_mid = veronoi_midline(origPoly.points,verbose)
    midpts = veronoi_mid.points.copy()
    midpts = midpts[np.argsort(midpts[:, 0])]
    if verbose:
        p = pv.Plotter()
        p.add_mesh(origPoly, color="black", opacity=0.1, label="origPoly")
        p.add_mesh(sortedPoly, color="orange", label="sortedPoly")
        p.add_mesh(veronoi_mid, color="green", label="veronoi_mid")
        p.add_legend()
        p.set_background("white")
        p.show()
    circles = []
    quads = []
    smashs = []
    edges = []
    farpts = []
    farptsids = []
    attempts = 0
    valid_checkPoints = []
    found_limits = {"low": False,
                    "high": False}
    for limit in found_limits.keys():

        while found_limits[limit] == False:
            attempts += 1

            if verbose:
                p = pv.Plotter()
                p.add_mesh(sortedPoly, color="orange", label="sortedPoly")
                #p.add_mesh(pv.PolyData(trypt), color="green", label="trypt")
                p.add_legend()
                p.set_background("white")
                p.show()

            while (found_limits[limit] != True):

                if limit == "low":
                    random_idx = np.random.randint(-int((0.15 + 0.15 * (attempts / 100)) * len(veronoi_mid.points)), -1)
                elif limit == "high":
                    random_idx = np.random.randint(0, int((0.15 + 0.15 * (attempts / 100)) * len(veronoi_mid.points)))

                trypt = midpts[random_idx]

                closest = closest_node_index(np.array([trypt[0], trypt[1], 0]), sortedPoly.points)
                closest_dist = vecAbs(np.array([trypt[0], trypt[1], 0]) - sortedPoly.points[closest])

                add = (attempts / 100) * closest_dist * np.array(
                    [2 * (-0.5 + np.random.rand()), 2 * (-0.5 + 1 * np.random.rand()), 0])

                shift_Try = trypt + add

                up = veronoi_mid.points[random_idx - 1][:]
                down = veronoi_mid.points[random_idx + 1][:]

                if up[0] > down[0]:
                    mid_tangent = (up - down)
                elif up[0] <= down[0]:
                    mid_tangent = (down - up)

                mid_angle = -angle_between(mid_tangent, np.array([0, 1, 0])) / np.pi * 180

                count_ang = 0
                while (found_limits[limit] != True and count_ang <= 10):
                    mid_angle += (attempts / 100) * np.random.randint(-15.5, 15.5)
                    count_ang += 1

                    try_center = np.array([shift_Try[0], shift_Try[1], 0])
                    try_radius = closest_dist + np.random.rand() * closest_dist * (0.15+0.15*(attempts / 100))
                    try_circle = pv.Cylinder(try_center,  # center
                                             (0, 0, 1),  # direction
                                             try_radius,  # radius
                                             closest_dist,  # height
                                             1000,  # resolution
                                             )

                    smash = linePoly.slice_along_line(polyline_from_points(try_circle.slice(normal="z").points))

                    try_quad = pv.Plane(center=try_center,
                                        i_size=2 * try_radius,
                                        j_size=2 * try_radius,
                                        direction=(0, 0, 1),
                                        i_resolution=1,
                                        j_resolution=1)
                    try_quad = try_quad.extract_feature_edges()
                    try_quad.translate(-try_center)
                    try_quad.rotate_z(mid_angle)
                    try_quad.translate(try_center)

                    if len(smash.points) == 4:

                        if verbose:
                            p = pv.Plotter()
                            p.add_mesh(sortedPoly, color="orange", label="sortedPoly")
                            p.add_mesh(pv.PolyData(trypt), color="green", label="trypt")
                            p.add_mesh(try_circle.slice(normal="z"), color="black", label="try_circle")
                            p.add_mesh(try_quad, color="black", label="try_quad")
                            p.add_mesh(smash, color="red", label="smash")
                            p.add_legend()
                            p.set_background("white")
                            p.show()

                        first_edge = pv.Line(try_quad.points[1], try_quad.points[2], 100)
                        second_edge = pv.Line(try_quad.points[3], try_quad.points[0], 100)
                        first_no = pv.Line(try_quad.points[0], try_quad.points[1], 100)
                        second_no = pv.Line(try_quad.points[2], try_quad.points[3], 100)

                        counter = np.zeros(4)

                        for smashpt in smash.points:
                            first_edge_dist = vecAbs(
                                first_edge.points[closest_node_index(smashpt, first_edge.points)] - smashpt)
                            second_edge_dist = vecAbs(
                                second_edge.points[closest_node_index(smashpt, second_edge.points)] - smashpt)
                            first_no_dist = vecAbs(
                                first_no.points[closest_node_index(smashpt, first_no.points)] - smashpt)
                            second_no_dist = vecAbs(
                                second_no.points[closest_node_index(smashpt, second_no.points)] - smashpt)

                            counter[
                                np.argmin([first_edge_dist, second_edge_dist, first_no_dist, second_no_dist])] += 1

                        if counter[0] == 2 and counter[1] == 2:

                            smashquad = smash.delaunay_2d()
                            smashquad = smashquad.extract_feature_edges()
                            smashquad_dist_first_edge = [
                                vecAbs(first_edge.points[closest_node_index(i, first_edge.points)] - i) for i in
                                smashquad.points]
                            smashquad_dist_second_edge = [
                                vecAbs(second_edge.points[closest_node_index(i, second_edge.points)] - i) for i in
                                smashquad.points]
                            smash_firstpts = []
                            smash_secpts = []
                            for idx, pt in enumerate(smashquad.points):
                                if smashquad_dist_first_edge[idx] > smashquad_dist_second_edge[idx]:
                                    smash_firstpts.append(pt)
                                else:
                                    smash_secpts.append(pt)
                            crosslines = [pv.Line(smash_firstpts[0], smash_secpts[0], 10),
                                          pv.Line(smash_firstpts[0], smash_secpts[1], 10),
                                          pv.Line(smash_firstpts[1], smash_secpts[0], 10),
                                          pv.Line(smash_firstpts[1], smash_secpts[1], 10)]

                            deletelines = []

                            for idx, cl in enumerate(crosslines):
                                if idx < 2:
                                    rg = (2, 4)
                                elif idx >= 2:
                                    rg = (0, 2)
                                for test_cl_idx in range(*rg):
                                    test_cl = crosslines[test_cl_idx]
                                    tester = pv.Line(test_cl.points[1], test_cl.points[-2])
                                    cross = tester.slice_along_line(cl)
                                    if cross.number_of_points > 0:
                                        deletelines.append(idx)

                            crosslines = [i for j, i in enumerate(crosslines) if j not in deletelines]
                            otherlines = [pv.Line(crosslines[0].points[0], crosslines[1].points[0], 2),
                                          pv.Line(crosslines[0].points[1], crosslines[1].points[1], 2)]

                            mids = [crosslines[0].points[5], crosslines[1].points[5]]

                            checkPoints = extract_edge_poi(try_center, try_radius, mids, limit, origPoly, verbose)

                            if len(checkPoints) == 0:
                                break
                            edges.append(first_no)
                            edges.append(second_no)
                            edges.append(first_edge)
                            edges.append(second_edge)

                            if any([i.length < try_radius / 2 for i in otherlines]):
                                break

                            farpt = [distant_node_index(i, checkPoints) for i in mids]

                            if verbose:
                                p = pv.Plotter()
                                p.add_mesh(sortedPoly, color="orange", label="sortedPoly")
                                p.add_mesh(pv.PolyData(trypt), color="green", label="trypt")
                                p.add_mesh(try_circle.slice(normal="z"), color="black", label="try_circle")
                                p.add_mesh(try_quad, color="black", label="try_quad")
                                p.add_mesh(smash, color="red", label="smash")
                                p.add_mesh(np.array([checkPoints[i] for i in farpt]), color="yellow",
                                           point_size=15, label="farpt")
                                p.add_mesh(np.array(mids), color="black", label="mids")
                                p.add_mesh(veronoi_mid, color="yellow", label="veronoi_mid")
                                p.add_legend()
                                p.set_background("white")
                                p.show()

                            if all_equal(farpt):

                                farpts.append(checkPoints[farpt[-1]])

                                oids = [np.where((sortedPoly.points == i).all(axis=1))[0][0] for i in farpts]

                                farptsids.append(oids[-1])
                                found_limits[limit] = True
                                smashs.append(smash)
                                quads.append(try_quad)
                                circles.append(try_circle)
                                valid_checkPoints.append(checkPoints)

                                if verbose:
                                    p = pv.Plotter()
                                    p.add_mesh(sortedPoly, color="orange", label="sortedPoly")
                                    p.add_mesh(pv.PolyData(trypt), color="green", label="trypt")
                                    p.add_mesh(try_circle.slice(normal="z"), color="black", label="try_circle")
                                    p.add_mesh(try_quad, color="black", label="try_quad")
                                    p.add_mesh(smash, color="red", label="smash")
                                    p.add_mesh(np.array([checkPoints[i] for i in farpt]), color="yellow",
                                               point_size=15, label="farpt")
                                   # p.add_mesh(realfarpt, color="red", point_size=20)
                                    p.add_mesh(np.array(mids), color="black", label="mids")
                                    p.add_mesh(veronoi_mid, color="yellow", label="veronoi_mid")
                                    #p.add_mesh(points[0],color="black",point_size=10)
                                    p.add_legend()
                                    p.set_background("white")
                                    p.show()

    ind_vk = farptsids[[i[0] for i in farpts].index(min([i[0] for i in farpts]))]
    ind_hk = farptsids[[i[0] for i in farpts].index(max([i[0] for i in farpts]))]
    return ind_hk, ind_vk, veronoi_mid


def extractSidePolys(ind_hk, ind_vk, sortedPoly):
    xs,ys = list(sortedPoly.points[::,0]), list(sortedPoly.points[::,1])

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
    psPoly = pv.PolyData(np.stack((x_ps, y_ps, np.zeros(len(x_ps)))).T)
    ssPoly = pv.PolyData(np.stack((x_ss, y_ss, np.zeros(len(x_ss)))).T)
    return psPoly, ssPoly


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
    xmids, ymids = midsPoly.points[::,0], midsPoly.points[::,1]
    vk_tangent = np.stack((xmids[0] - xmids[1], ymids[0] - ymids[1], 0)).T
    hk_tangent = np.stack((xmids[-2] - xmids[-1], ymids[-2] - ymids[-1], 0)).T
    camber_angle_vk = angle_between(vk_tangent, np.array([0, 1, 0])) / np.pi * 180
    camber_angle_hk = angle_between(hk_tangent, np.array([0, 1, 0])) / np.pi * 180
    return camber_angle_hk, camber_angle_vk


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


def getBoundaryValues(x_bounds, y_bounds):
    x = x_bounds
    y = y_bounds

    x_min = min(x)
    x_max = max(x)

    x_inlet = []
    x_outlet = []

    y_inlet = []
    y_outlet = []
    x_peri = []
    y_peri = []

    def sort_value(y, u):

        y = np.asarray(y)
        u = np.asarray(u)

        idx = np.argsort(y)

        y = np.asarray(y)[idx]
        u = np.asarray(u)[idx]

        new_y = []
        new_u = []

        for i in range(len(y)):
            new_y.append(y[i])
            new_u.append(u[i])

        return new_u

    def sort_value2(y, u):

        y = np.asarray(y)
        u = np.asarray(u)

        idx = np.argsort(y)

        y = np.asarray(y)[idx]
        u = np.asarray(u)[idx]

        new_y = []
        new_u = []

        for i in range(len(y)):
            new_y.append(y[i])
            new_u.append(u[i])

        return new_y, new_u

    for i in range(len(x)):
        if x[i] < x_min + 0.0000001 and x[i] > x_min - 0.0000001:
            x_inlet.append(x[i])
            y_inlet.append(y[i])
        if x[i] < x_max + 0.0000001 and x[i] > x_max - 0.0000001:
            x_outlet.append(x[i])
            y_outlet.append(y[i])
        if x[i] < x_max and x[i] > x_min:
            x_peri.append(x[i])
            y_peri.append(y[i])

    y_inlet, x_inlet = sort_value2(y_inlet, x_inlet)

    y_outlet, x_outlet = sort_value2(y_outlet, x_outlet)

    x_max_upper = x_outlet[-1]
    x_max_lower = x_outlet[0]

    x_peri.append(x_outlet[-1])
    y_peri.append(y_outlet[-1])
    x_peri.append(x_outlet[0])
    y_peri.append(y_outlet[0])

    x_upper_peri = []
    y_upper_peri = []

    x_upper_peri.append(x_inlet[-1])
    y_upper_peri.append(y_inlet[-1])

    x_lower_peri = []
    y_lower_peri = []

    x_lower_peri.append(x_inlet[0])
    y_lower_peri.append(y_inlet[0])

    x_start = -99999

    while x_start != x_max_upper:

        dists = []

        for i in range(len(x_peri)):
            dist = np.sqrt((x_upper_peri[-1] - x_peri[i]) ** 2 + (y_upper_peri[-1] - y_peri[i]) ** 2)
            dists.append(dist)

        x_peri = sort_value(dists, x_peri)
        y_peri = sort_value(dists, y_peri)

        x_upper_peri.append(x_peri[0])
        y_upper_peri.append(y_peri[0])

        x_start = x_upper_peri[-1]

        del x_peri[0]
        del y_peri[0]

    x_upper_peri.append(x_outlet[-1])
    y_upper_peri.append(y_outlet[-1])

    x_start = -99999

    while x_start != x_max_lower:

        dists = []

        for i in range(len(x_peri)):
            dist = np.sqrt((x_lower_peri[-1] - x_peri[i]) ** 2 + (y_lower_peri[-1] - y_peri[i]) ** 2)
            dists.append(dist)

        x_peri = sort_value(dists, x_peri)
        y_peri = sort_value(dists, y_peri)

        x_lower_peri.append(x_peri[0])
        y_lower_peri.append(y_peri[0])

        x_start = x_lower_peri[-1]

        del x_peri[0]
        del y_peri[0]

    x_lower_peri.append(x_outlet[0])
    y_lower_peri.append(y_outlet[0])

    return y_inlet, x_inlet, y_outlet, x_outlet, x_lower_peri, y_lower_peri, x_upper_peri, y_upper_peri


def getGeom2DVTUSLice2(mesh, alpha):
    midspan_slice, midspan_z = slice_midspan_z(mesh)
    midspan_slice = midspan_slice.compute_normals()
    geo = midspan_slice.extract_feature_edges()

    # points_complete = alle punkte auf dem mittelschnitt mit domain
    points_complete = midspan_slice.points

    points_bounds = geo.points

    x_outer_bounds, y_outer_bounds = calcConcaveHull(points_complete[:, 0], points_complete[:, 1], alpha)
    points_outer_bounds = np.stack((np.array(x_outer_bounds), np.array(y_outer_bounds)), axis=-1)

    x_profil = []
    y_profil = []


    for i in range(len(points_bounds)):
        if np.array([points_bounds[i][0], points_bounds[i][1]]) not in points_outer_bounds:
            x_profil.append(points_bounds[i][0])
            y_profil.append(points_bounds[i][1])

    return x_outer_bounds, y_outer_bounds, x_profil, y_profil, midspan_z


def rotatePoints(origin, x, y, angle):
    """
    Rotate a point counterclockwise by a given angle around a given origin.
    """

    angle = angle * (math.pi / 180.0)

    ox, oy = origin

    new_x = []
    new_y = []

    for i in range(len(x)):
        qx = ox + math.cos(angle) * (x[i] - ox) - math.sin(angle) * (y[i] - oy)
        qy = oy + math.sin(angle) * (x[i] - ox) + math.cos(angle) * (y[i] - oy)

        new_x.append(qx)
        new_y.append(qy)

    return new_x, new_y


def GetProfileValuesMidspan(case):
    midspan_slice, midspan_z = case.get_midspan_z()
    midspan_slice = midspan_slice.compute_normals()
    geo = midspan_slice.extract_feature_edges()

    # points_complete = alle punkte auf dem mittelschnitt mit domain
    points_complete = midspan_slice.points

    points_bounds = geo.points
    x_outer_bounds, y_outer_bounds = calcConcaveHull(points_complete[:, 0], points_complete[:, 1],
                                                     alpha=case.CascadeCoeffs.alpha)

    points_outer_bounds = np.stack((np.array(x_outer_bounds), np.array(y_outer_bounds)), axis=-1)

    x_profil = []
    y_profil = []

    indexes_profil_points = []

    for i in range(len(points_bounds)):
        if np.array([points_bounds[i][0], points_bounds[i][1]]) not in points_outer_bounds:
            indexes_profil_points.append(i)
            x_profil.append(points_bounds[i][0])
            y_profil.append(points_bounds[i][1])

    x_ss, y_ss, x_ps, y_ps = sortProfilePoints(x_profil, y_profil, alpha=case.CascadeCoeffs.alpha)

    indexes_ss = []
    indexes_ps = []

    x_l_ax_ss = []
    x_l_ax_ps = []

    for i in range(len(x_ss)):
        x_l_ax_ss.append((x_ss[i] - min(x_ss)) / (max(x_ss) - min(x_ss)))
        new_index = np.where(points_bounds == [x_ss[i], y_ss[i], 0])[0][0]
        indexes_ss.append(new_index)

    for i in range(len(x_ps)):
        x_l_ax_ps.append((x_ps[i] - min(x_ps)) / (max(x_ps) - min(x_ps)))
        new_index = np.where(points_bounds == [x_ps[i], y_ps[i], 0])[0][0]
        indexes_ps.append(new_index)

    assert len(
        indexes_ps) > 0, "sortProfilePoints did not detect a pressureside. Maybe adjust case.CascadeCoeffs.alpha ? "
    assert len(
        indexes_ss) > 0, "sortProfilePoints did not detect a suctionside. Maybe adjust case.CascadeCoeffs.alpha ? "

    values_ss = []
    values_ps = []
    value_names = []

    wall_shear_stress_ss = []
    wall_shear_stress_ps = []
    wall_shear_stress_explike_ss = []
    wall_shear_stress_explike_ps = []

    normals_ps = np.asarray([geo.cell_normals[i] for i in indexes_ps])

    for i in geo.point_arrays:
        value_names.append(i)
        values_ss.append([])
        values_ps.append([])
        for idx in indexes_ps:
            values_ps[-1].append(geo.point_arrays[i][idx])
        for idx in indexes_ss:
            values_ss[-1].append(geo.point_arrays[i][idx])

    for i in range(len(indexes_ss)):

        dudx = values_ss[value_names.index('dudx')][i]
        dudy = values_ss[value_names.index('dudy')][i]
        dudz = values_ss[value_names.index('dudz')][i]
        dvdx = values_ss[value_names.index('dvdx')][i]
        dvdy = values_ss[value_names.index('dvdy')][i]
        dvdz = values_ss[value_names.index('dvdz')][i]
        dwdx = values_ss[value_names.index('dwdx')][i]
        dwdy = values_ss[value_names.index('dwdy')][i]
        dwdz = values_ss[value_names.index('dwdz')][i]

        if i > 0:

            face_normal_delta_x = x_ss[i] - x_ss[i - 1]
            face_normal_delta_y = y_ss[i] - y_ss[i - 1]

        else:

            face_normal_delta_x = x_ss[i + 1] - x_ss[i]
            face_normal_delta_y = y_ss[i + 1] - y_ss[i]

        [face_normal_delta_x], [face_normal_delta_y] = rotatePoints([0, 0], [face_normal_delta_x],
                                                                    [face_normal_delta_y], -90.0)
        face_normal = np.array([face_normal_delta_x, face_normal_delta_y, 0])
        face_normal = face_normal / np.linalg.norm(face_normal)
        rho = values_ss[value_names.index('rho')][i]

        if 'nu' not in value_names:
            nu = Sutherland_Law(values_ss[value_names.index('T')][i])
        else:
            nu = values_ss[value_names.index('nu')][i]
        p = values_ss[value_names.index('p')][i]

        wall_shear_stress_vec = calcWallShearStress(dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, face_normal,
                                                    rho, nu, p)
        wall_shear_stress_abs = np.linalg.norm(wall_shear_stress_vec)
        if wall_shear_stress_vec[0] < 0 and face_normal[0] > 0:  # hier ist noch ein fehler.
            # Auch der normalen Vektor des faces muss mit einbezogen werden
            wall_shear_stress_abs = -wall_shear_stress_abs
        elif wall_shear_stress_vec[0] > 0 and face_normal[0] < 0:
            wall_shear_stress_abs = -wall_shear_stress_abs

        wall_shear_stress_explike_ss.append(np.sqrt(wall_shear_stress_vec[0] ** 2 + wall_shear_stress_vec[1] ** 2))
        wall_shear_stress_ss.append(wall_shear_stress_abs)

    for i in range(len(indexes_ps)):
        dudx = values_ps[value_names.index('dudx')][i]
        dudy = values_ps[value_names.index('dudy')][i]
        dudz = values_ps[value_names.index('dudz')][i]
        dvdx = values_ps[value_names.index('dvdx')][i]
        dvdy = values_ps[value_names.index('dvdy')][i]
        dvdz = values_ps[value_names.index('dvdz')][i]
        dwdx = values_ps[value_names.index('dwdx')][i]
        dwdy = values_ps[value_names.index('dwdy')][i]
        dwdz = values_ps[value_names.index('dwdz')][i]
        face_normal = -normals_ps[i]
        rho = values_ps[value_names.index('rho')][i]

        if 'nu' not in value_names:
            nu = Sutherland_Law(values_ps[value_names.index('T')][i])
        else:
            nu = values_ps[value_names.index('nu')][i]
        p = values_ps[value_names.index('p')][i]

        wall_shear_stress_vec = calcWallShearStress(dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, face_normal,
                                                    rho, nu, p)

        wall_shear_stress_abs = np.linalg.norm(wall_shear_stress_vec)
        if wall_shear_stress_vec[0] < 0 and face_normal[
            0] > 0:  # hier ist noch ein fehler. Auch der normalen Vektor des faces muss mit einbezogen werden
            wall_shear_stress_abs = -wall_shear_stress_abs
        elif wall_shear_stress_vec[0] > 0 and face_normal[0] < 0:
            wall_shear_stress_abs = -wall_shear_stress_abs

        # wall_shear_stress_ss.append(wall_shear_stress_abs)

        wall_shear_stress_ps.append(wall_shear_stress_abs)
        wall_shear_stress_explike_ps.append(np.sqrt(wall_shear_stress_vec[0] ** 2 + wall_shear_stress_vec[1] ** 2))

    values_ss.insert(0, x_l_ax_ss)
    values_ps.insert(0, x_l_ax_ps)
    values_ss.insert(0, y_ss)
    values_ps.insert(0, y_ps)
    values_ss.insert(0, x_ss)
    values_ps.insert(0, x_ps)
    values_ss.append(wall_shear_stress_ss)
    values_ss.append(wall_shear_stress_explike_ss)
    values_ps.append(wall_shear_stress_ps)
    values_ps.append(wall_shear_stress_explike_ps)

    value_names.insert(0, "x<sub>Ax</sub> / l<sub>Ax</sub>")
    value_names.insert(0, 'Y')
    value_names.insert(0, 'X')
    value_names.append('wall shear stress')
    value_names.append('wall shear stress exp like')

    return [value_names, [values_ss, values_ps]]


def getPitchValuesB2BSliceComplete(case, x):
    mesh = case.mesh_loaded_dict["fluid"]

    cut_plane = mesh.slice(normal="x", origin=(x, 0, 0))

    points = cut_plane.points

    npts = len(points)
    xx = np.zeros(npts)
    y = np.zeros(npts)
    zz = np.zeros(npts)

    for i in range(npts):
        pt = points[i]
        y[i] = pt[1]
        xx[i] = pt[0]
        zz[i] = pt[2]

    array_names = cut_plane.point_arrays.keys()
    values = []

    for arrname in array_names:
        array_values = cut_plane.point_arrays[arrname]

        y2, [array_values] = sort_values_by_pitch(y, [array_values])
        values.append(array_values)

    return y2, array_names, values


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


def closest_node_index(node, nodes):
    closest_index = distance.cdist([node], nodes).argmin()
    return closest_index


def distant_node_index(node, nodes):
    closest_index = distance.cdist([node], nodes).argmax()
    return closest_index


def inside_poly(polygon, points):
    """
    :param polygon: (x,y)
    :param points: (x,y)
    :return: list of True or False for indicating weather inside or not
    """
    path = mpltPath.Path(polygon)
    inside = path.contains_points(points)
    return inside


def all_equal(iterable):
    return iterable.count(iterable[0]) == len(iterable)


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
