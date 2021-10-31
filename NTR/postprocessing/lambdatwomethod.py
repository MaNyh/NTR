import pyvista as pv
import numpy as np

from NTR.utils.pyvista_utils import load_mesh
mesh = load_mesh(r"D:\NTR\examples\ChannelCase_les\03_Solution\VTK\02_Simcase_36000.vtk")

gradient = mesh.compute_derivative("U","qcriterion")
#gradient.rename_array("gradient","dU")
gradient.set_active_scalars("qcriterion")
gradient = gradient.contour(10)
a = gradient.threshold_percent([0.4, 0.55])
a.plot()
