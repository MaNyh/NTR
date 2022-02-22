import numpy as np
import pyvista as pv
import os
import vtk

from NTR.database.case_dirstructure import casedirs
from NTR.utils.mesh_handling.vtk_utils import cgnsReader

def translate_meshnames(mesh):
    array_names = mesh.array_names
    for name in array_names:
        if name == "Velocity":
            mesh.rename_array("Velocity", "U")
        if name == "Pressure":
            mesh.rename_array("Pressure", "p")
        if name == "Density":
            mesh.rename_array("Density", "rho")
        if name == "Temperature":
            mesh.rename_array("Temperature", "T")
    mesh = mesh.cell_data_to_point_data()
    return mesh

def load_mesh(path_to_mesh):
    assert os.path.isfile(path_to_mesh), path_to_mesh + " is not a valid file"
    extension = os.path.splitext(path_to_mesh)[1]
    if extension == ".vtk" or ".vtu":
        try:
            mesh = pv.UnstructuredGrid(path_to_mesh)
        except:
            try:
                mesh = pv.PolyData(path_to_mesh)
            except:
                print("error loading ", path_to_mesh)

    elif extension == ".cgns":
        cgns = cgnsReader(path_to_mesh)
        if str(type(cgns)).find('vtkStructuredGrid') != -1:
            mesh = pv.StructuredGrid(cgns)
        elif str(type(cgns)).find('vtkUnstructuredGrid') != -1:
            mesh = pv.UnstructuredGrid(cgns)
        elif str(type(cgns)).find('vtkMultiBlockDataSet') != -1:

            appendFilter = vtk.vtkAppendFilter()
            # Points with same coordinates are merged
            # with tolerance 0.0000001 GMC GLOBAL Properties
            appendFilter.MergePointsOn()
            appendFilter.SetTolerance(0.0000001)

#            mesh = pv.UnstructuredGrid()
            multiBlockMesh = pv.MultiBlock(cgns)
            for domainId in range(multiBlockMesh.GetNumberOfBlocks()):
                domain = multiBlockMesh.GetBlock(domainId)
                for blockId in range(domain.GetNumberOfBlocks()):
                    block = domain.GetBlock(blockId)
                    appendFilter.AddInputData(block)
                    appendFilter.Update()

            vtkmesh = appendFilter.GetOutput()
            mesh = pv.UnstructuredGrid(vtkmesh)

    elif extension == ".vtm":
        mesh = pv.PolyData()
        multiBlockMesh = pv.MultiBlock(path_to_mesh)
        for domainId in range(multiBlockMesh.GetNumberOfBlocks()):
            domain = multiBlockMesh.GetBlock(domainId)
            for blockId in range(domain.GetNumberOfBlocks()):
                block = domain.GetBlock(blockId)
                for patchId in range(block.GetNumberOfBlocks()):
                    patch = block.GetBlock(patchId)
                    mesh = mesh.merge(patch)

    mesh = translate_meshnames(mesh)
    return mesh


def mesh_scalar_gradients(mesh, array_name):
    gradients_mesh = mesh.compute_derivative(scalars=array_name)
    keys = np.array(["dudx", "dudy", "dudz", "dvdx", "dvdy", "dvdz", "dwdx", "dwdy", "dwdz"])
    keys = keys.reshape((3, 3))[:, :gradients_mesh["gradient"].shape[1]].ravel()
    gradients = dict(zip(keys, gradients_mesh["gradient"].T))
    for k, v in gradients.items():
        mesh[k] = v

    return mesh


def slice_midspan_z(mesh):
    bounds = mesh.bounds
    midspan_z = (bounds[5] - bounds[4]) / 2
    slice = mesh.slice(normal="z", origin=(0, 0, midspan_z))
    return slice, midspan_z


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
    cells = np.full((len(points) - 1, 3), 2, dtype=np.int_)
    cells[:, 1] = np.arange(0, len(points) - 1, dtype=np.int_)
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


def calc_dist_from_surface(surface_primary, surface_secondary, verbose=False):
    """
    Distance Between Two Surfaces / A Surface and a Pointcloud
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Compute the average thickness between two surfaces.

    We can compute the thickness between the two surfaces using a few different
    methods. We will demo a method where we compute the normals of the
    bottom surface, and then project a ray to the top surface to compute the
    distance along the surface normals.
    :param surface_primary: pv.UnstructuredGrid
    :param surface_secondary: pv.UnstructuredGrid / pv.PolyData
    :param verbose: plots?
    :return: surface_primary with distances from secondary
    """

    ###############################################################################
    # Ray Tracing Distance
    # ++++++++++++++++++++
    #

    h0n = surface_primary.compute_normals(point_normals=True, cell_normals=False,
                                          auto_orient_normals=True)

    ###############################################################################
    # Travel along normals to the other surface and compute the thickness on each
    # vector.

    h0n["distances"] = np.empty(surface_primary.n_points)
    for i in range(h0n.n_points):
        p = h0n.points[i]
        vec = h0n["Normals"][i] * h0n.length
        p0 = p - vec
        p1 = p + vec
        ip, ic = surface_secondary.ray_trace(p0, p1, first_point=True)
        dist = np.sqrt(np.sum((ip - p) ** 2))
        h0n["distances"][i] = dist

    if any(h0n["distances"] == 0):
        # Replace zeros with nans
        mask = h0n["distances"] == 0
        h0n["distances"][mask] = np.nan
        np.nanmean(h0n["distances"])

    if verbose:
        ###############################################################################
        p = pv.Plotter()
        p.add_mesh(h0n, scalars="distances", smooth_shading=True)
        p.add_mesh(surface_secondary, color=True, opacity=0.75, smooth_shading=True)
        p.show()

    return h0n


def plot_geometry_tofile(path_to_sim, probes_to_plot, geometry_plots, plotname, zoom=1, point_size=6):
    pv.set_plot_theme("document")

    p = pv.Plotter(off_screen=True, window_size=[4800, 4800])
    probe_colors = ["red", "blue", "green", "yellow"]
    for probename, probepoly in probes_to_plot.items():
        p.add_mesh(probepoly, label=probename, point_size=point_size, color=probe_colors.pop(0))

    for geomname, geompoly in geometry_plots.items():
        p.add_mesh(geompoly)

    p.add_legend(bcolor=(1, 1, 1), )
    p.camera.position = (0, 0, 1)
    p.camera.roll += 270
    p.camera.zoom(zoom)
    p.show(screenshot=os.path.join(path_to_sim, "..", casedirs["data"], plotname))


def print_meshquality(mesh):
    mqual = mesh.compute_cell_quality()
    mquals = mqual["CellQuality"]
    print(mesh)
    print("mean: " + str(round(np.mean(mquals),2)))
    print("min: " + str(round(np.min(mquals),2)))
    print("max: " + str(round(np.max(mquals),2)))
