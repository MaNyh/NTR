import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

from NTR.utils.pyvista_utils import load_mesh


def vol_to_line(vtkmesh, array_name, ave_direction, verbose=False):
    """
    this function is assuming a structured grid without curved gridlines
    it extracts layers and averages them. currently the face-normals have to be alligned with the global coordinate system

    :param settings:
    :param verbose:
    :return:
    """
    mesh = vtkmesh
    dirs = {"x": 0, "y": 2, "z": 4}
    interpol_dir = dirs[ave_direction]

    rest = mesh.copy()
    pts = []
    meanvals = []
    pbar = tqdm(total=mesh.number_of_cells)

    while (rest.number_of_cells > 0):
        if verbose:
            p = pv.Plotter()
            p.add_mesh(mesh, opacity=0.5)
            p.add_mesh(rest)
            p.show()
        centers = rest.cell_centers()
        bounds = centers.bounds
        bnd = bounds[interpol_dir]

        ids = np.where(
            np.equal(centers.points[::, int(interpol_dir - interpol_dir / 2)], np.ones(len(centers.points)) * bnd))

        ids_negative = np.where(
            np.not_equal(centers.points[::, int(interpol_dir - interpol_dir / 2)], np.ones(len(centers.points)) * bnd))

        if len(ids_negative[0]) > 0:
            rest = rest.extract_cells(np.array([i for i in range(len(centers.points)) if not np.isin(i, ids[0])]))
        else:
            rest = pv.UniformGrid()
        layer = mesh.extract_cells(ids[0])
        mean = np.average(layer[array_name], axis=0)
        pts.append(bnd)
        meanvals.append(mean)
        pbar.update(layer.number_of_points)
    pbar.close()
    pos = np.array(pts)
    vals = np.array(meanvals)

    return pos, vals


vtkmesh = load_mesh(r"D:\CodingProjects\NTR\examples\ChannelCase_les\03_Solution\mittelung99Domain_148000.vtk")

array_name = "UPrime2Mean"
line_direction = "y"
pos, vals = vol_to_line(vtkmesh, array_name, line_direction, verbose=False)
pos, vals2 = vol_to_line(vtkmesh, "turbulenceProperties:RMean", line_direction, verbose=False)
vals += vals2
# vals können auch tensoren oder vektoren sein. hier nur für skalare getestet
plt.plot(pos, vals)
