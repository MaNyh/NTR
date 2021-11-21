import os

from NTR.utils.filehandling import yaml_dict_read
from NTR.utils.mesh_handling.pyvista_utils import load_mesh
settings_yml = "case_settings.yml"
settings = yaml_dict_read(settings_yml)
case_path = os.path.dirname(settings_yml)
path_to_volmesh = os.path.join(case_path, settings["post_settings"]["volmesh"])

mesh = load_mesh(path_to_volmesh)
surface = mesh.extract_surface()
regions = surface.connectivity()
