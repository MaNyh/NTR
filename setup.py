#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = ['Click>=7.0',
                "pip>=21.0.1",
                "bump2version>=0.5.11",
                "wheel>=0.33.6",
                "watchdog>=0.9.0",
                "flake8>=3.7.8",
                "tox>=3.14.0",
                "coverage>=4.5.4",
                "Sphinx>=1.8.5",
                "twine>=1.14.0",
                "matplotlib>=3.3.4",
                "numpy>=1.20.0",
                "vtk>=9.0.1",
                "scipy>=1.6.1",
                "pyvista>=0.25.0",
                "PyYaml>=5.4.1",
                "scikit-image>=0.18.3",
                "pandas>=1.3.3"
                ]

setup_requirements = [ ]

test_requirements = [ ]

setup(
    author="Numerical Test Rig",
    author_email='nyhuis@tfd.uni-hannover.de',
    python_requires='>=3.5',
    classifiers=[
        'Development Status :: 2 - Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.8',
    ],
    description="Numerical Test Rig. A Python-Package for Pre- and Postprocessing Computational Fluid Dynamics Simulations with a focus on transient and scale resolving simulations.",
    entry_points={
        'console_scripts': [
            'NTR=NTR.cli:main',
        ],
    },
    install_requires=requirements,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='NTR',
    name='NTR',
    packages=find_packages(include=['NTR', 'NTR.*']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/nyhuis/NTR',
    version='0.1.4',
    zip_safe=False,
)
