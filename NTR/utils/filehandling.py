import csv
import os
import pickle
import yaml
from functools import reduce


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


def write_pickle_protocolzero(file, args):
    with open(file, "wb") as Fobj:
        pickle.dump(args, Fobj, protocol=0)

def write_pickle(file, args):
    with open(file, "wb") as Fobj:
        pickle.dump(args, Fobj)

def read_pickle(file):
    with open(file, 'rb') as f:
        args = pickle.load(f)
    return args

def write_yaml_dict(fpath,data):
    with open(fpath, 'w') as outfile:
        yaml.dump(data, outfile, default_flow_style=False)


def readtxtfile(path_to_file):
    basepath = os.path.abspath(os.path.dirname(__file__))
    with open(os.path.join(basepath, os.path.basename(path_to_file)), "r") as fobj:
        content = fobj.readlines()
    return content

def write_file(path_to_file,content):
    with open(path_to_file, "w", newline='\n') as fobj:
        fobj.writelines(content)


def get_template_contents(templatepath, file_templates):
    template_contents = {}
    for key, vallist in file_templates.items():
        template_contents[key] = {}
        for val in vallist:
            template_contents[key][val] = "".join(readtxtfile(os.path.join(templatepath, key, val)))

    return template_contents

def get_directory_structure(rootdir):
    """
    Creates a nested dictionary that represents the folder structure of rootdir
    """
    dir = {}
    rootdir = rootdir.rstrip(os.sep)
    start = rootdir.rfind(os.sep) + 1
    for path, dirs, files in os.walk(rootdir):
        folders = path[start:].split(os.sep)
        subdir = dict.fromkeys(files)
        parent = reduce(dict.get, folders[:-1], dir)
        parent[folders[-1]] = subdir
    return dir

def walk_file_or_dir(root):
    if os.path.isfile(root):
        dirname, basename = os.path.split(root)
        yield dirname, [], [basename]
    else:
        for path, dirnames, filenames in os.walk(root):
            yield path, dirnames, filenames
