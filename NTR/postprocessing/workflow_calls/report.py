from NTR.utils.mesh_handling.pyvista_utils import load_mesh
import pyvista as pv
import numpy as np
import matplotlib.pylab as pl
from matplotlib.colors import ListedColormap


def create_transparent_cmap(cmap,low,high):
    """
    :param cmap: matplotlib.pylab.cm - colormap
    :return: cap - transparent
    """
    norm = max([abs(low),abs(high)])
    my_cmap = cmap(np.arange(cmap.N))
    # Set alpha
    my_cmap[:, -1] = abs(np.linspace(-low/norm, high/norm, cmap.N))
    # Create new colormap
    my_cmap = ListedColormap(my_cmap)
    return my_cmap


input = r"Z:\Arbeit\Promotion_TFD_2020\02_TuV\06_Simulationsdaten\RAS\02_SEC\Virginiatech_Sensivitaeten\01_Results\secdev-6400-12\output\cgns\TRACE.cgns"

mesh = load_mesh(input)
pv.set_plot_theme("document")


annotations = {
    1: "+100%",
    0: "Unchanged",
    -1: "-100%",
}

res = 4800
title_size = int(0.02 * res)
sargs = dict(
    title_font_size=title_size,
    label_font_size=int(0.016 * res),
    shadow=True,
    n_labels=3,
    italic=True,
    # fmt="%.1f",
    font_family="arial",
)
differences = mesh.copy()
differences.clear_data()


ref = input.replace("secdev", "reference").replace("12", "10")
refmesh = load_mesh(ref)


for an in ["TurbulentEnergyKinetic"]:
    differences[an] = (refmesh[an] - mesh[an]) / refmesh[an]


low = min(differences["TurbulentEnergyKinetic"])
high = max(differences["TurbulentEnergyKinetic"])
differences["TKE_normed"] = differences["TurbulentEnergyKinetic"]/max([high,abs(low)])

my_cmap_RdBu = create_transparent_cmap(pl.cm.seismic,low,high)
#my_cmap_binary = create_transparent_cmap(pl.cm.binary)



differencespts = pv.PolyData(differences.points)
differencespts["TKE_normed"]=differences.point_data["TKE_normed"]
#edges = differences.extract_all_edges()
feature_edges = mesh.extract_feature_edges()
p = pv.Plotter(off_screen=True)
p.add_mesh(differencespts,annotations=annotations, scalars="TKE_normed", cmap=my_cmap_RdBu , scalar_bar_args=sargs)
p.add_mesh(feature_edges,show_scalar_bar=False)
#p.add_volume(differences, cmap=my_cmap_RdBu, opacity="sigmoid")

p.show(window_size=[res, res], title=input)
