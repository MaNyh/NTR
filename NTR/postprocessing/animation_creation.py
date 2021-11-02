# -*- coding: utf-8 -*-
"""
Created on Mon May  3 00:33:31 2021

@author: malte
"""


import pyvista as pv
import imageio
import os
from tqdm import tqdm

from NTR.database.case_dirstructure import casedirs
from NTR.utils.filehandling import yaml_dict_read


def create(path_to_yaml_dict):
    settings = yaml_dict_read(path_to_yaml_dict)
    casepath = os.path.abspath(os.path.dirname(path_to_yaml_dict))

    cutplanepath = os.path.join(casepath,casedirs["solution"],"postProcessing","cuttingPlane")
    assert os.path.isdir(cutplanepath),"no data-directory found"
    assert "post_settings" in settings.keys(), "no post_settings defined"
    assert "animation_creation" in settings["post_settings"].keys(), "no animation_creation defined"
    assert "variable" in settings["post_settings"]["animation_creation"].keys(), "no variable for animation defined"
    assert "resolution_x" in settings["post_settings"]["animation_creation"].keys(), "no resolution for animation defined"
    assert "resolution_y" in settings["post_settings"]["animation_creation"].keys(), "no resolution for animation defined"
    assert "low_scale_limit" in settings["post_settings"]["animation_creation"].keys(), "no low_scale_limit for animation defined"
    assert "high_scale_limit" in settings["post_settings"]["animation_creation"].keys(), "no high_scale_limit for animation defined"
    var = settings["post_settings"]["animation_creation"]["variable"]
    resolution_x = int(settings["post_settings"]["animation_creation"]["resolution_x"])
    resolution_y = int(settings["post_settings"]["animation_creation"]["resolution_y"])

    cpos = settings["post_settings"]["animation_creation"]["cpos"]

    low_scale = settings["post_settings"]["animation_creation"]["low_scale_limit"]
    high_scale = settings["post_settings"]["animation_creation"]["high_scale_limit"]

    dirs = [i for i in os.listdir(cutplanepath) if os.path.isdir(os.path.join(cutplanepath,i))]
    vtkname = var+"_constantPlane.vtk"
    plotter = pv.Plotter(notebook=False, off_screen=True)

    plotter.open_gif(os.path.join(casepath,casedirs["data"],var + "_animation.gif"))

    for d in tqdm(dirs):
        target = os.path.join(cutplanepath,d,vtkname)
        plotter.theme = pv.themes.DocumentTheme()
        mesh = pv.PolyData(target)
        mesh.rotate_z(90)
        actor_mesh = plotter.add_mesh(mesh,cmap="coolwarm")
        plotter.update_scalar_bar_range((low_scale, high_scale))
        plotter.show_axes()

        plotter.render()
        plotter.write_frame()
        plotter.remove_actor(actor_mesh)
        plotter.clear()

    plotter.close()
