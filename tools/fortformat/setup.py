#!/usr/bin/env python3
from distutils.core import setup


setup(
    name='fortformat',
    version='0.2',
    description='Basic Fortnet Input Format Class',
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
Basic Fortnet Input Format Class
--------------------------------
This basic Python class implements the Fortnet input file format. The input
geometries, stored in a list of ASE atoms objects, and targets, stored in a
list of Numpy arrays, conveniently get dumped to disk as simple text files.
''',
    requires=['numpy']
)
