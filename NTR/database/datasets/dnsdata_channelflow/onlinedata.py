import numpy as np
import requests
import tempfile
import os
import h5py
from scipy import integrate
import matplotlib.pyplot as plt

import NTR
from NTR.utils.filehandling import yaml_dict_read

# todo: check where to get these values online
nu = 3.50000e-04  # DNS-Moser
u_tau = 6.37309e-02  # DNS-Moser


def get(case_name):
    basedir = os.path.join(os.path.dirname(NTR.__file__), "database", "datasets", "dnsdata_channelflow")
    cases_urls = yaml_dict_read(os.path.join(basedir, "cases.yml"))
    case_links = cases_urls[case_name]
    tempdir = tempfile.TemporaryDirectory()
    data = {}
    for var, file in case_links.items():
        extension = os.path.splitext(file)[1]

        tmpfile = os.path.join(tempdir.name, var)
        r = requests.get(file, allow_redirects=True)
        print("loaded " + file)
        if extension == ".dat":
            with open(tmpfile, "w") as varf:
                varf.write(r.content.decode("ascii"))
            data[var] = np.genfromtxt(tmpfile, comments='%')
        elif extension == ".h5":
            with open(tmpfile, "wb") as varf:
                varf.write(r.content)
            with h5py.File(tmpfile, "r") as hdf:
                # List all groups
                h5pydat = {}
                ls = list(hdf.keys())
                for item in ls:
                    h5pydat[item] = np.array(hdf.get(item))
            data[var] = h5pydat
    return data


def create_bcdata_from_dnsdatabase(case_name,verbose=False):
    simdata = get(case_name)
    DNS_mean = simdata["vel_pressure_mean_prof"]
    Euu = simdata["energy_spec"]["Euu_kx"]

    ###################################### Returning spatial-two-point-correlation & length scale from given spectra
    dx_range = np.linspace(0, 4.0 * np.pi, len(Euu))  # delta x space for two-point-correlation

    Ls = []
    Rs = [simdata["uu_prof"],simdata["uv_prof"],simdata["vv_prof"],simdata["uu_prof"]*0,simdata["ww_prof"]]
    Us = simdata["vel_pressure_mean_prof"][::,2]*u_tau
    Ps = simdata["vel_pressure_mean_prof"][::,5]*u_tau
    Ruu_returned = []
    for Euu_pt in Euu.T:
        Ruu_return_tmp = np.fft.ifft(((Euu_pt)) / 2)
        Ruu_return_real_tmp = (Ruu_return_tmp.real)
        Ruu_returned.append(Ruu_return_real_tmp)
        L_return_tmp = integrate.simps(Ruu_return_real_tmp / Ruu_return_real_tmp[0],
                                       dx_range)  # normalizing with Ruu_return_real_tmp[0]
        Ls.append(L_return_tmp)
    Ls = np.asarray(Ls)

    if verbose:
        plt.figure()
        plt.plot(Ls, DNS_mean[:, 0], '-o', label=case_name + ' - DNS-Moser 2015')
        plt.legend()
        plt.xlabel("L")
        plt.ylabel("y")
        plt.grid()

    return Ls, Rs, Us, Ps


create_bcdata_from_dnsdatabase("ReTau5200")
