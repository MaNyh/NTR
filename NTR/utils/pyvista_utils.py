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
    return slice , midspan_z


def polyline_from_points(points):
    poly = pv.PolyData()
    poly.points = points
    the_cell = np.arange(0, len(points), dtype=np.int_)
    the_cell = np.insert(the_cell, 0, len(points))
    poly.lines = the_cell
    return poly


def lines_from_points(points):
    """Given an array of points, make a line set"""
    poly = pv.PolyData()
    poly.points = points
    cells = np.full((len(points)-1, 3), 2, dtype=np.int_)
    cells[:, 1] = np.arange(0, len(points)-1, dtype=np.int_)
    cells[:, 2] = np.arange(1, len(points), dtype=np.int_)
    poly.lines = cells
    return poly
