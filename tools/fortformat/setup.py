#!/usr/bin/env python3
from distutils.core import setup


setup(
    name='fortformat',
    version='0.4',
    description='Basic Fortnet IO Format Classes',
    author='T. W. van der Heide',
    url='https://github.com/vanderhe/fortnet',
    platforms='platform independent',
    package_dir={'': 'src'},
    packages=[''],
    classifiers=[
        'Programming Language :: Python',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering',
    ],
    long_description='''
Basic Fortnet IO Format Classes
-------------------------------
These basic Python classes implement the Fortnet input and output file format.
The Fnetdata class enables to create compatible HDF5 datasets, whereas the
Fnetout class extracts certain properties of the HDF5 output for later analysis.
''',
    requires=['numpy', 'h5py']
)
