import csv
import os
import pickle

import yaml


def yaml_dict_read(yml_file):
    args_from_yaml = {}

    with open(yml_file, "r", newline='') as Fobj:
        document = yaml.load_all(Fobj, Loader=yaml.FullLoader)
        for settings in document:
            for key, value in settings.items():
                args_from_yaml[key] = value
    return args_from_yaml


def read_csv(csv_filepath):
    with open(csv_filepath,"r", newline='') as csvobj:
        spamreader = csv.reader(csvobj, delimiter='\t', quotechar='|')
        data = []
        for row in spamreader:
            data.append(row)
    return data


def write_pickle(file, args):
    with open(file, "wb") as Fobj:
        pickle.dump(args, Fobj, protocol=0)


def write_yaml_dict(fpath,data):
    with open(fpath, 'w') as outfile:
        yaml.dump(data, outfile, default_flow_style=False)


def readtxtfile(path_to_file):
    basepath = os.path.abspath(os.path.dirname(__file__))
    with open(os.path.join(basepath,path_to_file), "r") as fobj:
        content = fobj.readlines()
    return content
