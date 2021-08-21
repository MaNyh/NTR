import numpy as np
from scipy.spatial import distance
from scipy.spatial.distance import squareform, pdist


def calc_largedistant_idx(x_koords, y_koords):
    A = np.dstack((x_koords, y_koords))[0]
    D = squareform(pdist(A))
    #    N = np.max(D)
    I = np.argmax(D)
    I_row, I_col = np.unravel_index(I, D.shape)

    index_1 = I_row
    index_2 = I_col

    return index_1, index_2


def closest_node_index(node, nodes):
    closest_index = distance.cdist([node], nodes).argmin()
    return closest_index


def distant_node_index(node, nodes):
    closest_index = distance.cdist([node], nodes).argmax()
    return closest_index
