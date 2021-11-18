import pyvista as pv
import numpy as np

from NTR.utils.mesh_handling.pyvista_utils import load_mesh
#mesh = load_mesh(r"C:\Users\Nyhuis\Desktop\Neuer Ordner\VTK\01_data_12000.vtk")

mesh = load_mesh(r"C:\Users\Nyhuis\Desktop\Neuer Ordner\VTK\mittelung99Domain_148000.vtk")
gradient = mesh.compute_derivative("U", "qcriterion")
gradient.set_active_scalars("qcriterion")
plotstuff = gradient.contour(8)
plotstuff = plotstuff.smooth(400)
#a = gradient.threshold_percent([0.49, 0.5])
#a.plot()
