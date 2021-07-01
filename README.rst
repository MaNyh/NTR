===
NTR
===


.. image:: https://img.shields.io/pypi/v/NTR.svg
        :target: https://pypi.python.org/pypi/NTR

.. image:: https://img.shields.io/travis/nyhuis/NTR.svg
        :target: https://travis-ci.com/nyhuis/NTR

.. image:: https://readthedocs.org/projects/NTR/badge/?version=latest
        :target: https://NTR.readthedocs.io/en/latest/?version=latest
        :alt: Documentation Status




Numerical Test Rig. A Python-Package for Pre- and Postprocessing Computational Fluid Dynamics Simulations with a focus on transient and scale resolving simulations.


* Free software: MIT license
* Documentation: https://NTR.readthedocs.io.


Installation
-------------

Install package using pip
'pip install PATH/TO/PACKAGE'

Or install via
'python setup.py install'

As a developer, install development-requirements via
'pip install -r requirements_dev.txt'


Be aware that this package is in early development. It is not refactored well and functions are not validated or tested correctly

Features
-------------
postprocessing
    -openfoam
    .createProfileData
preprocessing
    -openfoam
        -create_probes

utils
    .aeroFunctions
    .boundaryLayerFunctions
    .case
    .casereader
    .coefficients
    .create_geom
    .pyvista_utils
    .simFunctions
    .solver_variable_dicts
    .thermoFunctions

Usage igg_cascademeshing
---------------------------

See ./examples/"meshing-case"
Create a case-directory

./case
    -pointcloudfile.txt

    -case_settings.yml

    -call_stuff.py


This structure is not final and must be optimized in the future. Set up settings.yml for your specific case. Right now, also the path to igg has to be set here


Usage meshing
---------------------------

Create meshes using the igg-python-interpreter. Settings from the "case_settings.yml" are transmitted using a pickle-object.



Usage create_case
---------------------------

Create cfd-cases with this function. Use an example-case to learn how to set the "case.yml".
The function will create standardized cases for conducting and monitoring simulations.
Only OpenFOAM-CascadeCases implemented yet. Large-Eddy-Simulations are running but not validated.
RAS-cases are not parametrized yet. The RAS-Package is only a dummy-package right now

Based on the settings, the function is calling create_probes (see below)




Usage create_probes
---------------------------

Create probe-dictionaries with this algorithm. You will need the profile-surface as a vtk (PolyData).
Set your options in the "case.yml"-file. Try out an example, you need a \*.vtk - domain for a case.
Only cascade-cases implemented.

The function can be called using

from NTR.preprocessing.openfoam import createProbesProfileDict

createProbesProfileDict("pathtovtk",pden_Probes_Profile_SS, pden_Probes_Profile_PS,
                            interval_time_steps_probes, output_path, tolerance=1e-6

The dictionary will be written out in the working directory

Usage testing-module
---------------------------

The module can be partially tested using pytest and the testing-module in /tests/tesT_NTR.py
There will be plenty of warnings, because things are not coded correctly.
The testing will perform tests on the main geometrical functions.
Further tests are welcome. If things are not working correctly, let me know and I will develop tests.


Credits
-------

The NTR-Package is based on the NTR-Package by Mark Ziesse (ziesse@tfd.uni-hannover.de)


This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
