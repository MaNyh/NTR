import numpy as np
import pyvista as pv
import os

def load_mesh(path_to_mesh):
    assert os.path.isfile(path_to_mesh),path_to_mesh + " is not a valid file"
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


def constructWallMesh(meshList):
    print("constructing surfacemesh from wall meshes ...")
    wallMesh = pv.UnstructuredGrid()

    for mesh in meshList:
        m = load_mesh(mesh)
        m = m.compute_normals()
        wallMesh = wallMesh.merge(m)

    return wallMesh


def calc_dist_from_surface(surface_one, surface_two, verbose=False):
    """
    Distance Between Two Surfaces / A Surface and a Pointcloud
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Compute the average thickness between two surfaces.

    We can compute the thickness between the two surfaces using a few different
    methods. We will demo a method where we compute the normals of the
    bottom surface, and then project a ray to the top surface to compute the
    distance along the surface normals.
    :param surface_one: pv.UnstructuredGrid
    :param surface_two: pv.UnstructuredGrid / pv.PolyData
    :param verbose: plots?
    :return: surface with distance
    """

    ###############################################################################
    # Ray Tracing Distance
    # ++++++++++++++++++++
    #

    h0n = surface_one.compute_normals(point_normals=True, cell_normals=False,
                                      auto_orient_normals=True)

    ###############################################################################
    # Travel along normals to the other surface and compute the thickness on each
    # vector.

    h0n["distances"] = np.empty(surface_one.n_points)
    for i in range(h0n.n_points):
        p = h0n.points[i]
        vec = h0n["Normals"][i] * h0n.length
        p0 = p - vec
        p1 = p + vec
        ip, ic = surface_two.ray_trace(p0, p1, first_point=True)
        dist = np.sqrt(np.sum((ip - p) ** 2))
        h0n["distances"][i] = dist

    if any(h0n["distances"]==0):
        # Replace zeros with nans
        mask = h0n["distances"] == 0
        h0n["distances"][mask] = np.nan
        np.nanmean(h0n["distances"])

    if verbose:
        ###############################################################################
        p = pv.Plotter()
        p.add_mesh(h0n, scalars="distances", smooth_shading=True)
        p.add_mesh(surface_two, color=True, opacity=0.75, smooth_shading=True)
        p.show()


    return surface_one
