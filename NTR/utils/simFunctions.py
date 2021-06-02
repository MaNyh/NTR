import numpy as np


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


def sort_value3(y, u):
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

    return new_u  #
