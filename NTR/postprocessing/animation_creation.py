# -*- coding: utf-8 -*-
"""
Created on Mon May  3 00:33:31 2021

@author: malte
"""


import pyvista as pv
import imageio
import os

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
    var = settings["post_settings"]["animation_creation"]["variable"]

    dirs = [i for i in os.listdir(cutplanepath) if os.path.isdir(os.path.join(cutplanepath,i))]
    vtkname = var+"_constantPlane.vtk"
    frame = "frame.jpg"

    with imageio.get_writer("animation.gif", mode='I') as writer:
        for d in dirs:
            target = os.path.join(cutplanepath,d,vtkname)
            plotter = pv.Plotter(off_screen=True)
            plotter.background_color = (1,1,1)

            mesh = pv.PolyData(target)
            mesh.rotate_z(90)
            #mesh.plot(cpos=[0,0,1])
            plotter.add_mesh(mesh,cmap="coolwarm")
            plotter.show_axes()
            plotter.update_scalar_bar_range((0,120))
            plotter.show(screenshot=frame, window_size=[800, 800],cpos=[0,0,1])
            plotter.close()
            image = imageio.imread(frame)
            writer.append_data(image)
            os.remove(frame)
        writer.close()

