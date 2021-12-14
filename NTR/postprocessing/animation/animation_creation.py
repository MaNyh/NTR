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
        focus = animationsettings["focus"]
        view = animationsettings["view"]
        if "postprocess" in animationsettings.keys():
            postprocess = animationsettings["post"]
        else:
            postprocess = None

        if "autoscale" in animationsettings.keys():
            autoscale = animationsettings["autoscale"]
        else:
            autoscale = None

        low_scale = animationsettings["low_scale_limit"]
        high_scale = animationsettings["high_scale_limit"]

        dirs = [i for i in os.listdir(cutplanepath) if os.path.isdir(os.path.join(cutplanepath, i))]
        vtkname = var + "_constantPlane.vtk"
        plotter = pv.Plotter(notebook=False, off_screen=True, window_size=([resolution_x, resolution_y]))

        plotter.open_gif(os.path.join(casepath, casedirs["data"], var + "_" + animationname + ".gif"))
        plotter.background_color = (1, 1, 1)

        if autoscale:
            dir = dirs[0]
            target = os.path.join(cutplanepath, dir, vtkname)
            mesh = pv.PolyData(target)
            if postprocess == "divergence":
                mesh = mesh.compute_derivative(scalars=var, divergence=True)
                low_scale = min(mesh["divergence"])#animationsettings["low_scale_limit"]
                high_scale = max(mesh["divergence"])#animationsettings["high_scale_limit"]

        for d in tqdm(dirs):
            target = os.path.join(cutplanepath, d, vtkname)
            # plotter.theme = pv.themes.DocumentTheme()
            mesh = pv.PolyData(target)
            mesh.set_active_scalars(var)
            if postprocess == "divergence":
                mesh = mesh.compute_derivative(scalars=var, divergence=True)

                mesh.set_active_scalars("divergence")
            actor_mesh = plotter.add_mesh(mesh, cmap="coolwarm")
            plotter.update_scalar_bar_range((low_scale, high_scale))
            plotter.camera_position = [cpos, focus, view]




            plotter.show_axes()

            plotter.render()
            # plotter.show(cpos=cpos)
            plotter.write_frame()
            plotter.remove_actor(actor_mesh)
            plotter.clear()

        plotter.close()
