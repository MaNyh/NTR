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


def sort_values_by_pitch(y, list_values):
    y = np.asarray(y)

    new_list_values = []

    for i in range(len(list_values)):
        list_values[i] = np.asarray(list_values[i])
        new_list_values.append([])

    idx = np.argsort(y)

    y = np.asarray(y)[idx]

    for i in range(len(list_values)):
        list_values[i] = np.asarray(list_values[i])[idx]

    new_y = []

    for i in range(len(y)):

        new_y.append(y[i])

        for j in range(len(list_values)):
            new_list_values[j].append(list_values[j][i])

    return new_y, new_list_values
