import math
import numpy as np
import pyvista as pv

from scipy.interpolate import UnivariateSpline
from scipy.spatial import Delaunay

from NTR.utils.pyvista_utils import mesh_scalar_gradients, slice_midspan_z
from NTR.utils.thermoFunctions import Sutherland_Law
from NTR.utils.boundaryLayerFunctions import calcWallShearStress
from NTR.utils.simFunctions import sort_values_by_pitch

def calcMidPoints(x1, y1, x2, y2):
    x_mid_ss = []
    y_mid_ss = []

    for i in range(len(x1)):

        dists = []

        for j in range(len(x2)):
            dist = ((x1[i] - x2[j]) ** 2 + (y1[i] - y2[j]) ** 2) ** (0.5)

            dists.append(dist)

        index_p = np.argmin(dists)

        p_x = x2[index_p]
        p_y = y2[index_p]

        def midpoint(x1, y1, x2, y2):
            return ((x1 + x2) / 2, (y1 + y2) / 2)

        x_mid, y_mid = midpoint(p_x, p_y, x1[i], y1[i])

        x_mid_ss.append(x_mid)
        y_mid_ss.append(y_mid)

    return x_mid_ss, y_mid_ss


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
    spl = UnivariateSpline(x_mcl, y_mcl, k=5)
    deri_spl = spl.derivative()

    deri_values = []

    deri_vk = np.sin(np.deg2rad(beta1 - 90.0)) / np.cos(np.deg2rad(beta1 - 90.0))
    deri_hk = -np.cos(np.deg2rad(beta2)) / np.sin(np.deg2rad(beta2))

    deri_diff_vk = []
    deri_diff_hk = []

    for i in range(len(x_mcl)):
        deri_values.append(deri_spl(x_mcl[i]))
        deri_diff_vk.append(((deri_vk - deri_values[-1]) ** 2) ** 0.5)
        deri_diff_hk.append(((deri_hk - deri_values[-1]) ** 2) ** 0.5)

    index_deri_vk = np.argmin(deri_diff_vk)
    index_deri_hk = np.argmin(deri_diff_hk)

    x_values = []
    y_values = []

    for i in range(len(x_mcl)):

        if index_deri_vk <= i <= index_deri_hk:
            x_values.append(x_mcl[i])
            y_values.append(y_mcl[i])

    delta_x_vk = x_values[0] - x_inlet
    delta_y_vk = np.tan(np.deg2rad(beta1 - 90.0)) * delta_x_vk

    p_inlet_x = x_values[0] - delta_x_vk
    p_inlet_y = y_values[0] - delta_y_vk

    delta_x_hk = x_outlet - x_values[-1]
    delta_y_hk = -delta_x_hk / np.tan(np.deg2rad(beta2))

    p_outlet_x = x_values[-1] + delta_x_hk
    p_outlet_y = y_values[-1] + delta_y_hk

    x_mpsl = [p_inlet_x] + x_values + [p_outlet_x]
    y_mpsl = [p_inlet_y] + y_values + [p_outlet_y]

    for i in range(len(x_mpsl)):
        y_mpsl[i] = y_mpsl[i] + 0.5 * t

    return x_mpsl, y_mpsl


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


def sortProfilePoints(x, y, alpha=0.007):
    x, y = calcConcaveHull(x, y, alpha=alpha)

    ind_vk = x.index(min(x))
    ind_hk = x.index(max(x))

    x_ss = x[ind_hk:ind_vk]
    y_ss = y[ind_hk:ind_vk]

    y_ps = y[:ind_hk-1] + y[ind_vk-1:]
    x_ps = x[:ind_hk-1] + x[ind_vk-1:]


    x_ps, y_ps = zip(*sorted(zip(x_ps, y_ps)))

    #x_ss, y_ss = zip(*sorted(zip(x_ss, y_ss)))
    return x_ss, y_ss, x_ps, y_ps


def rotate_points(origin, x, y, angle):
    """
    Rotate a point counterclockwise (rotation around z-axis) by a given angle around a given origin.
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


def calcMeanCamberLine(x, y, beta1, beta2):
    # vk und hk bestimmen

    x, y = zip(*sorted(zip(x, y)))

    ind_vk, ind_hk = calc_vk_hk(x, y, beta1, beta2)

    x_vk = x[ind_vk]
    y_vk = y[ind_vk]

    x_hk = x[ind_hk]
    y_hk = y[ind_hk]

    x_ss, y_ss, x_ps, y_ps = sortPoints(x, y, ind_vk, ind_hk)

    x_mid_ss = []
    y_mid_ss = []

    x_mid_ss, y_mid_ss = calcMidPoints(x_ss, y_ss, x_ps, y_ps)
    x_mid_ps, y_mid_ps = calcMidPoints(x_ps, y_ps, x_ss, y_ss)
    x_mids, y_mids = calcMidPoints(x_mid_ps, y_mid_ps, x_mid_ss, y_mid_ss)

    return x_mids, y_mids, x_ss, y_ss, x_ps, y_ps, x_vk, y_vk, x_hk, y_hk


def sortPoints(x, y, ind_vk, ind_hk):
    # Punkte der Saugseite bestimmen

    # ersten Punkt nach der Vorderkante bestimmen

    dists = []
    indizes = []

    for i in range(len(x)):
        if i != ind_vk and y[i] > y[ind_vk]:
            dists.append(((x[i] - x[ind_vk]) ** 2 + (y[i] - y[ind_vk]) ** 2) * (0.5))
            indizes.append(i)

    ind_ss_p2 = indizes[dists.index(min(dists))]

    indizes = list(range(len(x)))

    indizes.remove(ind_vk)
    indizes.remove(ind_ss_p2)

    indizes_ss = []

    indizes_ss.append(ind_vk)
    indizes_ss.append(ind_ss_p2)

    ind = ind_ss_p2

    while ind != ind_hk:

        dists = []
        inds = []

        point = (x[ind], y[ind])
        for i in range(len(indizes)):
            point2 = (x[indizes[i]], y[indizes[i]])
            dist = ((point2[0] - point[0]) ** 2 + (point2[1] - point[1]) ** 2) * (0.5)

            if indizes[i] not in indizes_ss:
                dists.append(dist)
                inds.append(indizes[i])

        indizes_ss.append(inds[dists.index(min(dists))])
        indizes.remove(inds[dists.index(min(dists))])
        ind = inds[dists.index(min(dists))]

        indizes_ps = list(indizes)
        indizes_ps.insert(0, ind_vk)
        indizes_ps.append(ind_hk)

    x_ss = []
    y_ss = []

    for i in range(len(indizes_ss)):
        x_ss.append(x[indizes_ss[i]])
        y_ss.append(y[indizes_ss[i]])

    x_ps = []
    y_ps = []

    for i in range(len(indizes_ps)):
        x_ps.append(x[indizes_ps[i]])
        y_ps.append(y[indizes_ps[i]])

    return x_ss, y_ss, x_ps, y_ps


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


def getGeom2DVTUSLice2(case):
    #mesh = pv.UnstructuredGrid(path_to_mesh)
    #cut_plane_polydata = slice_midspan_z(mesh)
    polyData = case.get_midspan_z()#cut_plane_polydata
    bounds = case.mesh_loaded_dict["fluid"].bounds
    midspan_z = (bounds[5]-bounds[4])/2

    # Boundary Edges / Zellen extrahieren
    featureEdges = polyData.extract_feature_edges()
    points_complete = featureEdges.points
    points_bounds = np.array([featureEdges.extract_cells(i).bounds for i in range(len(featureEdges.points))])

    x_outer_bounds, y_outer_bounds = calcConcaveHull(points_complete[:, 0], points_complete[:, 1], case.CascadeCoeffs.alpha)
    points_outer_bounds = np.stack((np.array(x_outer_bounds), np.array(y_outer_bounds)), axis=-1)

    x_profil = []
    y_profil = []

    indexes_profil_points = []

    for i in range(len(points_bounds)):
        if np.array([points_bounds[i][0], points_bounds[i][1]]) not in points_outer_bounds:
            indexes_profil_points.append(i)
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

    midspan_slice = case.get_midspan_z()
    midspan_slice = midspan_slice.compute_normals()
    geo = midspan_slice.extract_feature_edges()

    #points_complete = alle punkte auf dem mittelschnitt mit domain
    points_complete = midspan_slice.points

    points_bounds = geo.points
    x_outer_bounds, y_outer_bounds = calcConcaveHull(points_complete[:, 0], points_complete[:, 1], alpha=case.CascadeCoeffs.alpha)

    points_outer_bounds = np.stack((np.array(x_outer_bounds), np.array(y_outer_bounds)), axis=-1)

    x_profil = []
    y_profil = []

    indexes_profil_points = []

    for i in range(len(points_bounds)):
        if np.array([points_bounds[i][0], points_bounds[i][1]]) not in points_outer_bounds:
            indexes_profil_points.append(i)
            x_profil.append(points_bounds[i][0])
            y_profil.append(points_bounds[i][1])

    #x_profil = points_outer_bounds[::,0]
    #y_profil = points_outer_bounds[::,1]
    profile_points = np.stack((np.array(x_profil), np.array(y_profil)), axis=-1)

    x_ss, y_ss, x_ps, y_ps = sortProfilePoints(x_profil, y_profil, alpha=case.CascadeCoeffs.alpha)

    indexes_ss = []
    indexes_ps = []

    x_l_ax_ss = []
    x_l_ax_ps = []

    for i in range(len(x_ss)):
        x_l_ax_ss.append((x_ss[i] - min(x_ss)) / (max(x_ss) - min(x_ss)))
        indexes_ss.append(indexes_profil_points[profile_points.tolist().index([x_ss[i], y_ss[i]])])

    for i in range(len(x_ps)):
        x_l_ax_ps.append((x_ps[i] - min(x_ps)) / (max(x_ps) - min(x_ps)))
        indexes_ps.append(indexes_profil_points[profile_points.tolist().index([x_ps[i], y_ps[i]])])

    assert len(indexes_ps) > 0, "sortProfilePoints did not detect a pressureside. Maybe adjust case.CascadeCoeffs.alpha ? "
    assert len(indexes_ss) > 0, "sortProfilePoints did not detect a suctionside. Maybe adjust case.CascadeCoeffs.alpha ? "


    values_ss = []
    values_ps = []
    value_names = []

    wall_shear_stress_ss = []
    wall_shear_stress_ps = []
    wall_shear_stress_explike_ss = []
    wall_shear_stress_explike_ps = []

    #normals_ss = np.asarray([midspan_slice.cell_normals[i] for i in indexes_ss])
    normals_ps = np.asarray([midspan_slice.cell_normals[i] for i in indexes_ps])

    cells_ps = midspan_slice.extract_cells(indexes_ps)
    cells_ss = midspan_slice.extract_cells(indexes_ss)

    for i in midspan_slice.array_names:
        values_ss.append(cells_ss[i])
        values_ps.append(cells_ps[i])
        value_names.append(i)

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
        rho = values_ss[value_names.index('rhoMean')][i]

        if 'nu' not in value_names:
            nu = Sutherland_Law(values_ss[value_names.index('TMean')][i])
        else:
            nu = values_ss[value_names.index('nuMean')][i]
        p = values_ss[value_names.index('pMean')][i]


        wall_shear_stress_vec = calcWallShearStress(dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, face_normal, rho, nu, p)
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
            nu = Sutherland_Law(values_ps[value_names.index('TMean')][i])
        else:
            nu = values_ps[value_names.index('nuMean')][i]
        p = values_ps[value_names.index('pMean')][i]

        wall_shear_stress_vec = calcWallShearStress(dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, face_normal,
                                                    rho, nu, p)
        wall_shear_stress_abs = np.linalg.norm(wall_shear_stress_vec)
        if wall_shear_stress_vec[0] < 0 and face_normal[0] > 0:  # hier ist noch ein fehler. Auch der normalen Vektor des faces muss mit einbezogen werden
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

    cut_plane = mesh.slice(normal="x",origin=(x,0,0))

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

    #noa = len(cut_plane.array_names)
    array_names = cut_plane.point_arrays.keys()
    values = []

    for arrname in array_names:


        array_values = cut_plane.point_arrays[arrname]


        y2, [array_values] = sort_values_by_pitch(y, [array_values])
        values.append(array_values)

        # daten nach y sortieren

    return y2, array_names, values
