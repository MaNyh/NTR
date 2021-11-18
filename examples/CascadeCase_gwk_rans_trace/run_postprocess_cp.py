from NTR.utils.mesh_handling.pyvista_utils import load_mesh

fname = r"D:\NTR\examples\CascadeCase_gwk_rans_trace\03_Solution\TRACE.cgns"

mesh = load_mesh(fname)
mesh.plot()
