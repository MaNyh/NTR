import math
import numpy as np

from scipy.optimize import minimize

from NTR.utils.geom_functions.pointcloud import calcConcaveHull
from NTR.utils.geom_functions.profileparas import sortProfilePoints
from NTR.utils.fluid_functions.thermoFunctions import Sutherland_Law
from NTR.utils.fluid_functions.boundaryLayerFunctions import calcWallShearStress
from NTR.utils.simFunctions import sort_values_by_pitch
from NTR.pyvista_utils import slice_midspan_z


def calcConcaveHull_optimize(xs, ys, alpha):
    def func(alpha):
        xss, yss = calcConcaveHull(xs, ys, alpha)
        return len(xs) - len(xss)

    res = minimize(func, alpha, method='Powell',
                   options={'xatol': 1e-8, 'disp': True})

    opt_x, opt_y = calcConcaveHull(xs, ys, res.x[0])
    return opt_x, opt_y


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
