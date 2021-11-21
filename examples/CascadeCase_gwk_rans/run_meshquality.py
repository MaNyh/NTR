from NTR.utils.mesh_handling.pyvista_utils import load_mesh, print_meshquality

path = r"D:\NTR\examples\CascadeCase_gwk_rans\03_Solution\VTK\02_Simcase_30000.vtk"
mesh = load_mesh(path)
print_meshquality(mesh)
