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
    frame = "frame.jpg"

    pv.set_plot_theme("document")
    with imageio.get_writer(os.path.join(casepath,casedirs["data"],var + "_animation.gif"), mode='I') as writer:
        for d in tqdm(dirs):
            target = os.path.join(cutplanepath,d,vtkname)
            plotter = pv.Plotter(off_screen=True)

            mesh = pv.PolyData(target)
            mesh.rotate_z(90)
            plotter.add_mesh(mesh,cmap="coolwarm")
            plotter.show_axes()
            plotter.update_scalar_bar_range((low_scale,high_scale))
            plotter.camera.position = cpos
            plotter.show(screenshot=frame, window_size=[resolution_x,resolution_y])
            plotter.close()
            image = imageio.imread(frame)
            writer.append_data(image)
            os.remove(frame)
        writer.close()

