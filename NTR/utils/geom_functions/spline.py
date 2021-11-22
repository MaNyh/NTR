import numpy as np
from scipy.interpolate import interp1d, UnivariateSpline

from NTR.utils.geom_functions.vector import midpoint


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
        if abs(dist_1 - dist_2) / (dist_1 + dist_2) < 0.05:  # , rtol=tolerance
            x_mid_ss.append(x_mid)
            y_mid_ss.append(y_mid)

    x_mid_ss, y_mid_ss = zip(*sorted(zip(x_mid_ss, y_mid_ss)))

    return np.array(x_mid_ss), np.array(y_mid_ss)


def splineCurvature(xx,yy):
    dx_dt = np.gradient(xx)
    dy_dt = np.gradient(yy)
    ds_dt = np.sqrt(dx_dt * dx_dt + dy_dt * dy_dt)
    velocity = np.array([[dx_dt[i], dy_dt[i]] for i in range(dx_dt.size)])

    tangent = np.array([1 / ds_dt] * 2).transpose() * velocity
    tangent_x = tangent[:, 0]
    tangent_y = tangent[:, 1]

    deriv_tangent_x = np.gradient(tangent_x)
    deriv_tangent_y = np.gradient(tangent_y)

    dT_dt = np.array([[deriv_tangent_x[i], deriv_tangent_y[i]] for i in range(deriv_tangent_x.size)])

    length_dT_dt = np.sqrt(deriv_tangent_x * deriv_tangent_x + deriv_tangent_y * deriv_tangent_y)

    normal = np.array([1 / length_dT_dt] * 2).transpose() * dT_dt
    d2s_dt2 = np.gradient(ds_dt)
    d2x_dt2 = np.gradient(dx_dt)
    d2y_dt2 = np.gradient(dy_dt)
    curvature = np.abs(d2x_dt2 * dy_dt - dx_dt * d2y_dt2) / (dx_dt * dx_dt + dy_dt * dy_dt) ** 1.5

    return curvature


def curvature_splines(x, y=None, error=0.0001):
    """Calculate the signed curvature of a 2D curve at each point
    using interpolating splines.
    Parameters
    ----------
    x,y: numpy.array(dtype=float) shape (n_points, )
         or
         y=None and
         x is a numpy.array(dtype=complex) shape (n_points, )
         In the second case the curve is represented as a np.array
         of complex numbers.
    error : float
        The admisible error when interpolating the splines
    Returns
    -------
    curvature: numpy.array shape (n_points, )
    Note: This is 2-3x slower (1.8 ms for 2000 points) than `curvature_gradient`
    but more accurate, especially at the borders.
    """

    # handle list of complex case
    if y is None:
        x, y = x.real, x.imag

    t = np.arange(x.shape[0])
    std = error * np.ones_like(x)

    fx = UnivariateSpline(t, x, k=4, w=1 / np.sqrt(std))
    fy = UnivariateSpline(t, y, k=4, w=1 / np.sqrt(std))

    x_ = fx.derivative(1)(t)
    x__ = fx.derivative(2)(t)
    y_ = fy.derivative(1)(t)
    y__ = fy.derivative(2)(t)
    curvature = (x_* y__ - y_* x__) / np.power(x_** 2 + y_** 2, 3 / 2)
    return curvature
