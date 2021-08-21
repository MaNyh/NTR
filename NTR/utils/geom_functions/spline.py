import numpy as np
from scipy.interpolate import interp1d

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
