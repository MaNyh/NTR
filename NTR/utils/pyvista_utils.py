import numpy as np


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


def absVec(vec):
    return (vec[0]**2+vec[1]**2+vec[2]**2)**0.5

def absvec_array(array):
    return [absVec(vec) for vec in array]
