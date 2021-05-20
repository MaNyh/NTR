import math
import numpy as np
from scipy.interpolate import UnivariateSpline
from scipy.spatial import Delaunay


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
    x_mcl, y_mcl = Touple
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


def calcConcaveHull(x, y, alpha=0.007):
    """
    DOCUMENTATION INCOMPLETE

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


def calcMeanCamberLine(x, y, beta1, beta2):
    # vk und hk bestimmen

    x, y = zip(*sorted(zip(x, y)))

    ind_vk, ind_hk = calc_vk_hk(x, y, beta1, beta2)

    x_vk = x[ind_vk]
    y_vk = y[ind_vk]

    x_hk = x[ind_hk]
    y_hk = y[ind_hk]

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

        indizes = range(len(x))

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

    x_ss, y_ss, x_ps, y_ps = sortPoints(x, y, ind_vk, ind_hk)

    x_mid_ss = []
    y_mid_ss = []

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


def getGeom2DVTUSLice2(path_midspan_slice):
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(path_midspan_slice)
    reader.Update()
    data_complete = reader.GetOutput()
    data_complete.BuildLinks()

    mapper = vtk.vtkPointDataToCellData()
    if vtk.VTK_MAJOR_VERSION <= 5:
        mapper.AddInput(data_complete)
    else:
        mapper.AddInputData(data_complete)
    mapper.Update()
    data_bounds = mapper.GetOutput()
    data_bounds.BuildLinks()

    surfaceFilter = vtk.vtkDataSetSurfaceFilter()
    if vtk.VTK_MAJOR_VERSION <= 5:
        surfaceFilter.SetInput(data_bounds);
    else:
        surfaceFilter.SetInputData(data_bounds);
    surfaceFilter.Update();
    data_bounds = surfaceFilter.GetOutput()
    data_bounds.BuildLinks()

    # Slice von unstrukturiert zu PolyData umwandeln
    appendFilter = vtk.vtkDataSetSurfaceFilter()
    if vtk.VTK_MAJOR_VERSION <= 5:
        appendFilter.SetInput(data_bounds)
    else:
        appendFilter.SetInputData(data_bounds)
    appendFilter.Update()

    data_bounds = appendFilter.GetOutput()
    data_bounds.BuildLinks()

    mapper = vtk.vtkCellDataToPointData()
    if vtk.VTK_MAJOR_VERSION <= 5:
        mapper.AddInput(data_bounds)
    else:
        mapper.AddInputData(data_bounds)
    mapper.Update()
    data_bounds = mapper.GetOutput()
    data_bounds.BuildLinks()

    polyData = vtk.vtkPolyData()
    polyData.ShallowCopy(data_bounds)

    # Boundary Edges / Zellen extrahieren
    featureEdges = vtk.vtkFeatureEdges()
    if vtk.VTK_MAJOR_VERSION <= 5:
        featureEdges.SetInput(polyData)
    else:
        featureEdges.SetInputData(polyData)
    featureEdges.BoundaryEdgesOn()
    featureEdges.FeatureEdgesOff()
    featureEdges.ManifoldEdgesOff()
    featureEdges.NonManifoldEdgesOff()
    featureEdges.Update()

    data_bounds = featureEdges.GetOutput()
    data_bounds.BuildLinks()

    points_complete = vtk_to_numpy(data_complete.GetPoints().GetData())
    points_bounds = vtk_to_numpy(data_bounds.GetPoints().GetData())

    x_outer_bounds, y_outer_bounds = calcConcaveHull(points_complete[:, 0], points_complete[:, 1])
    points_outer_bounds = np.stack((np.array(x_outer_bounds), np.array(y_outer_bounds)), axis=-1)

    x_profil = []
    y_profil = []

    indexes_profil_points = []

    for i in range(len(points_bounds)):
        if np.array([points_bounds[i][0], points_bounds[i][1]]) not in points_outer_bounds:
            # print('yes')
            indexes_profil_points.append(i)
            x_profil.append(points_bounds[i][0])
            y_profil.append(points_bounds[i][1])

    # profile_points=np.unique(np.stack((np.array(x_profil), np.array(y_profil)), axis=-1) , axis=0)
    # profile_points=np.stack((np.array(x_profil), np.array(y_profil)), axis=-1)

    return x_outer_bounds, y_outer_bounds, x_profil, y_profil
