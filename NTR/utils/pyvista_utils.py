import numpy as np
import pyvista as pv

def load_mesh(path_to_mesh):
    try:
        mesh = pv.UnstructuredGrid(path_to_mesh)
    except:
        try:
            mesh = pv.PolyData(path_to_mesh)
        except:
            print("error loading ", path_to_mesh)

    print(mesh)
    return mesh

def mesh_scalar_gradients(mesh,array_name):
    gradients_mesh = mesh.compute_derivative(scalars=array_name)
    keys = np.array(["dudx", "dudy", "dudz", "dvdx", "dvdy", "dvdz", "dwdx", "dwdy", "dwdz"])
    keys = keys.reshape((3, 3))[:, :gradients_mesh["gradient"].shape[1]].ravel()
    gradients = dict(zip(keys, gradients_mesh["gradient"].T))
    for k, v in gradients.items():
        mesh[k] = v

    return mesh

def slice_midspan_z(mesh):
    bounds = mesh.bounds
    midspan_z = (bounds[5]-bounds[4])/2
    slice = mesh.slice(normal="z",origin=(0,0,midspan_z))
    return slice


