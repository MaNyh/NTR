import pyvista as pv
import numpy as np

from NTR.utils.pyvista_utils import load_mesh

mesh = load_mesh(r"D:\NTR\examples\ChannelCase_les\03_Solution\VTK\02_Simcase_160.vtk")
gradient = mesh.compute_derivative("U", "qcriterion")
gradient.set_active_scalars("qcriterion")
plotstuff = gradient.contour(12)
a = plotstuff.threshold_percent([0.4, 0.55])
a.plot()
