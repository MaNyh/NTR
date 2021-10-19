import numpy as np
import pyvista as pv
from matplotlib import path as mpltPath
from scipy.spatial import Delaunay
from scipy.interpolate import griddata,UnivariateSpline
from skimage.morphology import skeletonize

from NTR.utils.geom_functions.spline import refine_spline
from NTR.utils.pyvista_utils import polyline_from_points


def skeletonize_poly(points, verbose=False):

    points2d = points[::, 0:2]
    res = 1000
    pointsclosed = np.array(list(points2d) + [points2d[-1]])

    origPoly = pv.PolyData(points)

    bounds = origPoly.bounds
    x = np.linspace(bounds[0], bounds[1], res)
    y = np.linspace(bounds[2], bounds[3], res)

    pts = []
    for xi in x:
        for yi in y:
            pts.append((xi, yi))
    pts = np.array(pts)

    insides = inside_poly(pointsclosed, pts).astype(int)
    pts3d = []
    for idx, ins in enumerate(insides):
        if ins == 1:
            pts3d.append([pts[::, 0][idx], pts[::, 1][idx], 0])

    xv, yv = np.meshgrid(x, y, sparse=False, indexing='ij')

    grid_z1 = griddata(pts, insides, (xv, yv), method="nearest")

    skeleton = skeletonize(grid_z1)

    px = []
    py = []
    for idr, r in enumerate(skeleton):
        for idz, z in enumerate(r):
            if z != False:
                px.append(x[idr])
                py.append(y[idz])

    x_mid_ss, y_mid_ss = zip(*sorted(zip(px, py)))

    x_mid_ss_ref, y_mid_ref = refine_spline(x_mid_ss, y_mid_ss, 4000)
    spl = UnivariateSpline( x_mid_ss_ref, y_mid_ref)
    xs = np.linspace(x_mid_ss_ref[0], x_mid_ss_ref[-1], 1000)
    xout,yout = xs, spl(xs)

    midline = np.stack((xout,yout, np.zeros(len(yout)))).T
    outspline = polyline_from_points(midline)

    if verbose:
        p = pv.Plotter()
        p.add_mesh(origPoly)
        p.add_mesh(points, color="black")
        p.add_mesh(outspline, color="yellow")
        p.show()
    return outspline


def inside_poly(polygon, points):
    """
    :param polygon: (x,y)
    :param points: (x,y)
    :return: list of True or False for indicating weather inside or not
    """
    path = mpltPath.Path(polygon)
    inside = path.contains_points(points)
    return inside


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

            A = (s * (s - a) * (s - b) * (s - c))
            if A > 0:
                area = np.sqrt(A)

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
