import numpy as np
import os

class profile_data:
    def __init__(self, name, type, points, ind_hk, ind_vk):
        self.name = name
        self.type = type
        self.points = points
        self.ind_hk = ind_hk
        self.ind_vk = ind_vk


naca0009profile = profile_data("naca0009", "symmetric", np.loadtxt(os.path.join(os.path.dirname(__file__),"naca0009.pts")), 0, 34)
