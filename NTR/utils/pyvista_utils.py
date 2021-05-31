import numpy as np


def mesh_scalar_gradients(mesh,array_name):
    gradients_mesh = mesh.compute_derivative(scalars=array_name)
    keys = np.array(["du/dx", "du/dy", "du/dz", "dv/dx", "dv/dy", "dv/dz", "dw/dx", "dw/dy", "dw/dz"])
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
