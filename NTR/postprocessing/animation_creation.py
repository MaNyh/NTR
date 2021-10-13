# -*- coding: utf-8 -*-
"""
Created on Mon May  3 00:33:31 2021

@author: malte
"""


import pyvista as pv
import imageio
import os

dirs = [i for i in os.listdir() if os.path.isdir(i)]
vtkname = "U_zCut.vtk"
frame = "frame.jpg"

with imageio.get_writer("animation.gif", mode='I') as writer:
    for d in dirs:
        plotter = pv.Plotter(off_screen=True)
        plotter.background_color = (1,1,1)
        mesh = pv.PolyData(os.path.join(d,vtkname))
        mesh.rotate_z(90)
    #    mesh.plot(cpos=[0,0,1])
        plotter.add_mesh(mesh,cmap="coolwarm")
        plotter.show_axes()
        plotter.update_scalar_bar_range((0,120))
        plotter.show(screenshot=frame, window_size=[800, 800],cpos=[0,0,1])
        plotter.close()
        image = imageio.imread(frame)
        writer.append_data(image)
        os.remove(frame)
    writer.close()

