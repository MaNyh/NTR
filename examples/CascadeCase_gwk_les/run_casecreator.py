from NTR.preprocessing.openfoam.create_cascadecase_les import create_cascadecase_les
from NTR.preprocessing.openfoam.create_probes import create_probe_dicts


create_cascadecase_les("preprocessing_settings.yml")
create_probe_dicts("probes_settings.yml")
