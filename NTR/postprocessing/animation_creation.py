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

    cutplanepath = os.path.join(casepath, casedirs["solution"], "postProcessing", "cuttingPlane")

    assert os.path.isdir(cutplanepath), "no data-directory found"
    assert "post_settings" in settings.keys(), "no post_settings defined"
    assert "animation_creation" in settings["post_settings"].keys(), "no animation_creation defined"

    for animationname, animationsettings in settings["post_settings"]["animation_creation"].items():

        assert "variable" in animationsettings.keys(), "no variable for animation defined"
        assert "resolution_x" in animationsettings.keys(), "no resolution for animation defined"
        assert "resolution_y" in animationsettings.keys(), "no resolution for animation defined"
        assert "low_scale_limit" in animationsettings.keys(), "no low_scale_limit for animation defined"
        assert "high_scale_limit" in animationsettings.keys(), "no high_scale_limit for animation defined"

        var = animationsettings["variable"]
        resolution_x = int(animationsettings["resolution_x"])
        resolution_y = int(animationsettings["resolution_y"])

        cpos = animationsettings["cpos"]

        low_scale = animationsettings["low_scale_limit"]
        high_scale = animationsettings["high_scale_limit"]

        dirs = [i for i in os.listdir(cutplanepath) if os.path.isdir(os.path.join(cutplanepath, i))]
        vtkname = var + "_constantPlane.vtk"
        frame = "frame.jpg"

        pv.set_plot_theme("document")
        with imageio.get_writer(os.path.join(casepath, casedirs["data"], var + "_animation.gif"), mode='I') as writer:
            for d in tqdm(dirs):
                target = os.path.join(cutplanepath, d, vtkname)
                plotter = pv.Plotter(off_screen=True)

                mesh = pv.PolyData(target)
                # mesh.rotate_z(90)
                plotter.add_mesh(mesh, cmap="coolwarm")
                plotter.show_axes()
                plotter.update_scalar_bar_range((low_scale, high_scale))
                plotter.camera.position = cpos
                plotter.camera.roll += 270
                plotter.show(screenshot=frame, window_size=[resolution_x, resolution_y])
                plotter.close()
                image = imageio.imread(frame)
                writer.append_data(image)
                os.remove(frame)
            writer.close()
