import numpy as np
import os

from NTR.utils.functions import yaml_dict_read

profile_types = ["symmetric", "lift"]
list_of_datasets = ["naca0009.yml"]

datasets = {}
for datasetname in list_of_datasets:
    datasetdict = yaml_dict_read(os.path.join(os.path.dirname(__file__),"profiles", datasetname))
    assert isinstance(datasetdict["name"],str), "no valid name found in yaml-dcit"
    name = datasetdict["name"]
    datasets[name] = {}

    assert datasetdict["profile_type"] in profile_types, "types allowed:" + str(profile_types)
    profile_type = datasetdict["profile_type"]
    datasets[name]["profile_type"] = profile_type

    pointsraw = datasetdict["points"]
    points_floats = [float(i) for i in pointsraw.split(",")]
    coordsraw = np.array(points_floats)
    pxx,pyy = coordsraw[0:][::2],coordsraw[1:][::2]
    points = np.stack((pxx,pyy)).T
    datasets[name]["points"] = points

    assert isinstance(datasetdict["ind_hk"], int), "ind_hk must be an integer"
    ind_hk = datasetdict["ind_hk"]
    datasets[name]["ind_hk"] = ind_hk

    assert isinstance(datasetdict["ind_vk"], int), "ind_vk must be an integer"
    ind_vk = datasetdict["ind_vk"]
    datasets[name]["ind_vk"] = ind_vk


