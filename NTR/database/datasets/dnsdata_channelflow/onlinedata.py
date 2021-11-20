import numpy as np
import requests
import tempfile
import os
import h5py

import NTR
from NTR.utils.filehandling import yaml_dict_read

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


simdata = get("ReTau0180")
