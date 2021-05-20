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
--------
Install package using pip
'pip install PATH/TO/PACKAGE'

Or install via
'python setup.py install'

As a developer, install development-requirements via
'pip install -r requirements_dev.txt'

Usage
--------

Meshing:
See ./examples/"meshing-case"
Create a case-directory

./case
    pointcloudfile.txt
    settings.yml
    run.py

Run run.py

This structure is not final and must be optimized in the future.

Features
--------

* TODO

Credits
-------

The NTR-Package is based on the NTR-Package by Mark Ziesse (ziesse@tfd.uni-hannover.de)


This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
