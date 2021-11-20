import numpy as np
import requests
import tempfile
import os
import h5py
from scipy import integrate

import NTR
from NTR.utils.filehandling import yaml_dict_read

#todo: check where to get these values online
nu = 3.50000e-04 # DNS-Moser
u_tau = 6.37309e-02 # DNS-Moser

def get(case_name):
    basedir = os.path.join(os.path.dirname(NTR.__file__),"database","datasets","dnsdata_channelflow")
    cases_urls = yaml_dict_read(os.path.join(basedir,"cases.yml"))
    case_links = cases_urls[case_name]
    tempdir = tempfile.TemporaryDirectory()
    data = {}
    for var, file in case_links.items():
        extension = os.path.splitext(file)[1]

        tmpfile = os.path.join(tempdir.name,var)
        r = requests.get(file, allow_redirects=True)
        print("loaded " + file)
        if extension==".dat":
            with open(tmpfile, "w") as varf:
                varf.write(r.content.decode("ascii"))
            data[var]=np.genfromtxt(tmpfile, comments='%')
        elif extension==".h5":
            with open(tmpfile,"wb") as varf:
                varf.write(r.content)
            with h5py.File(tmpfile, "r") as hdf:
                # List all groups
                h5pydat = {}
                ls = list(hdf.keys())
                for item in ls:
                    h5pydat[item] = np.array(hdf.get(item))
            data[var]=h5pydat
    return data


def create_bcdata_from_dnsdatabase(case_name):
    simdata = get(case_name)
    DNS_mean = simdata["vel_pressure_mean_prof"]
    DNS_vel_fluc = simdata["covar_vel_prof"]
    Euu = simdata["energy_spec"]["Euu_kx"]
    kx = simdata["energy_spec"]["kx"]
    # TKE is normalized by 1/utau^2
    TKE_norm = 0.5 * (DNS_vel_fluc[:, 2] + DNS_vel_fluc[:, 3] + DNS_vel_fluc[:, 4])
    TKE = TKE_norm * u_tau ** 2  # adding dimension via u_tau
#    DNS_UMean = DNS_mean[:, 2] * u_tau

    ###################################### Returning spatial-two-point-correlation & length scale from given spectra
    dx_range_180 = np.linspace(0, 4.0 * np.pi, len(Euu))  # delta x space for two-point-correlation

    L_returned = []
    Ruu_returned = []
    for Euu_pt in Euu.T:
        Ruu_return_tmp = np.fft.ifft(((Euu_pt)) / 2)
        Ruu_return_real_tmp = (Ruu_return_tmp.real)
        Ruu_returned.append(Ruu_return_real_tmp)
        L_return_tmp = integrate.simps(Ruu_return_real_tmp / Ruu_return_real_tmp[0],dx_range_180)  # normalizing with Ruu_return_real_tmp[0]
        L_returned.append(L_return_tmp)
    L_returned = np.asarray(L_returned)

    plt.figure()
    plt.plot(L_returned, DNS_mean[:, 0], '-o', label=case_name + ' - DNS-Moser 2015')
    plt.legend()
    plt.xlabel("Lx")
    plt.ylabel("y/delta")
    plt.grid()


create_bcdata_from_dnsdatabase("ReTau1000")


