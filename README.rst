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


Be aware that this package is in constant development.
A test-module should keep things running. But not everything is tested and not everything can be tested within NTR itself.
When using this package, please contribute to this collection of methods. Write issues, write your own code and make merge-requests for methods, that bring benefits


Features
-------------
preprocessing
    * geometrycreation
    * experimental implementation of a naca-profile-generator
    * meshing for cascade-cases (using numeca igg and it's strong smoothing algorithm)
    * casecreation from ascii-based templates, independend from the solver that will be used
postprocessing
    * tbd

Examples
---------------------------

Currently not all examples are working properly. The implementation of the examples in a test-module should ensure, that this is the case in the future

Most cases will work though. Try to find a case that is suitable for your work and try to figure out how you can make use of the NTR methods

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


Creating probes
---------------------------

there is a method implemented for the creation of probes in openfoam-srs-cases

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
