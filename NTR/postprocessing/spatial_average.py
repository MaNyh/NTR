import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt

values = np.linspace(0, 10, 1000).reshape((20, 5, 10))
values.shape

# Create the spatial reference
grid = pv.UniformGrid()

# Set the grid dimensions: shape + 1 because we want to inject our values on
#   the CELL data
grid.dimensions = np.array(values.shape) + 1

# Edit the spatial reference
grid.origin = (100, 33, 55.6)  # The bottom left corner of the data set
grid.spacing = (1, 5, 2)  # These are the cell sizes along each axis

# Add the data values to the cell data
#grid["values"] = values.flatten(order="F")  # Flatten the array!
grid["U"]=np.ones(grid.number_of_cells)

settings = {"mesh": grid,
            "dir_1": "y"}


def vol_to_line(settings, verbose=True):
    """
    this function is assuming a structured grid without curved gridlines
    it extracts layers and averages them. currently the face-normals have to be alligned with the global coordinate system

    :param settings:
    :param verbose:
    :return:
    """
    mesh = settings["mesh"]
    array_name = "U"
    dirs = {"x":0,"y":2,"z":4}
    interpol_dir = dirs[settings["dir_1"]]

    rest = mesh.copy()
    pts = []
    meanvals = []
    while(rest.number_of_cells>0):
        if verbose:
            p = pv.Plotter()
            p.add_mesh(mesh, opacity=0.5)
            p.add_mesh(rest)
            p.show()
        centers = rest.cell_centers()
        bounds = centers.bounds
        bnd = bounds[interpol_dir]
        ids = np.where(np.equal(centers.points[::,int(interpol_dir-interpol_dir/2)],np.ones(len(centers.points))*bnd))
        ids_negative = np.array([i for i in range(len(centers.points)) if not np.isin(i,ids)])
        if len(ids_negative)>0:
            rest = rest.extract_cells(np.array([i for i in range(len(centers.points)) if not np.isin(i,ids)]))
        else:
            rest = pv.UniformGrid()
        layer = mesh.extract_cells(ids)
        mean = layer[array_name]
        pts.append(bnd)
        meanvals.append(mean)
    pos = np.array(pts)
    vals = np.array(meanvals)

    return pos, vals


pos, vals = vol_to_line(settings)

# vals können auch tensoren oder vektoren sein. hier nur für skalare getestet
plt.plot(pos,vals)
