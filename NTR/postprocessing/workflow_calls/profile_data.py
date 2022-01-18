import os
import matplotlib.pyplot as plt

from NTR.utils.filehandling import write_pickle, read_pickle
from NTR.utils.mesh_handling.pyvista_utils import load_mesh
from NTR.postprocessing.turbo.createProfileData import createProfileData

def profile_data_workflow(meshpath,output):

    # =============================================================================
    # Daten Einlesen
    # =============================================================================
    #profile_data = {}

    print(meshpath)
    refmesh = load_mesh(meshpath)

    midspan_z = (refmesh.bounds[5] - refmesh.bounds[4]) / 2
    alpha = 0.005
    post_slice_1_x = refmesh.bounds[0] + 1e-6
    post_slice_2_x = refmesh.bounds[1] - 1e-6
    kappa = 1.4
    As = 1.458e-06
    Ts = 110.4
    R_L = 287
    outflow = refmesh.slice(normal="x", origin=(post_slice_2_x, 0, 0)).compute_cell_sizes()
    outflow = outflow.point_data_to_cell_data()
    p_k = outflow["p"] * outflow["Area"] / (sum(outflow["Area"]))
    cp = 1004.5
    l = 0.13
    output_path = os.path.dirname(meshpath)

    profile_data = createProfileData(refmesh, midspan_z, alpha, post_slice_1_x, post_slice_2_x,
                                                 output_path,
                                                 kappa, R_L,
                                                 p_k, As, l, cp, Ts)

    """if not os.path.isdir(os.path.join(basedir,"02_PostProcessing")):
        os.mkdir(os.path.join(basedir,"02_PostProcessing"))
    """
    write_pickle(output,profile_data)


def plot_profiledata(basedir,resultpath, cascade_solution):
    fig, axs = plt.subplots(1, 3)

    for respath in cascade_solution:
        picklpath = os.path.join(basedir,resultpath,respath+"_profile_data.pkl")
        res = read_pickle(picklpath)
   #     print(res["mred"])
        axs[0].plot(res["mred"], res["inte_p2"]/res["inte_p1"],marker = "x")
        axs[1].plot(res["mred"], res["eta_is"],marker = "x" )
        axs[2].plot(res["mred"], res["delta_beta"],marker = "x")
    axs[0].set_xlabel('mred')
    axs[0].set_ylabel('Pi')
    axs[1].set_xlabel('mred')
    axs[1].set_ylabel('eta_is')
    axs[2].set_xlabel('mred')
    axs[2].set_ylabel('lift')
    #plt.show()
    plt.savefig("03_Plots/parastud.pdf")
    plt.close()

"""
def plot_profiledata(basedir,resultpath, cascade_solution):
    fig, axs = plt.subplots(1, 1)

    if not os.path.isdir("03_Plots"):
        os.mkdir("03_Plots")
        print("created dir")

    for respath in cascade_solution:
        picklpath = os.path.join(basedir,resultpath,respath+"_profile_data.pkl")
        res = read_pickle(picklpath)
        print(res["mred"])
        axs[0].plot(res["mred"], res["inte_p2"]/res["inte_p1"],marker = "x")
            axs[0].set_xlabel('mred')
    plt.savefig("03_Plots/profile_pressure.pdf")
    plt.close()
"""
