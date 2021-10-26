.. _sec-fortformat:
.. highlight:: none

##########
Fortformat
##########

[Input: `recipes/fortformat/`]

This chapter should serve as a tutorial guiding you through your first dataset
creation and data extraction after running Fortnet by using two Python classes,
provided by the `Fortformat <https://github.com/vanderhe/fortnet-python>`_
Python package:

.. toctree::
   :maxdepth: 1

   fnetdata.rst
   fnetout.rst

After these tutorials, you will be able to create a Fortnet compatible dataset,
based on the output files of your simulation package of choice (e.g. VASP),
as well as extracting the resulting output, i.e. predictions of the network.

*******************************************
Fortformat: Basic Fortnet IO Format Classes
*******************************************

This Python package provides two basic classes that implement the associated
Fortnet input and output file format. The input features and targets, stored in
lists of Numpy arrays, conveniently get dumped to disk as HDF5 files, while
Fortnet's output may be extracted from the HDF5 output files as written by a
validation or prediction run.


Installation
============

Please note, that this package has been tested for Python 3.X support. Its usage
additionally requires

  - `numerical Python <https://numpy.org/doc/stable/reference/>`_ (`numpy`)
  - `pythonic HDF5 <http://www.h5py.org/>`_ (`h5py`)
  - `Atomic Simulation Environment <https://wiki.fysik.dtu.dk/ase/>`_ (`ase`)

as well as the `pytest` framework in order to run the regression tests.

Via the Python Package Index
----------------------------

First, make sure you have an up-to-date version of pip installed::

  python -m pip install --upgrade pip

The package can be downloaded and installed via pip into the active Python
interpreter (preferably using a virtual python environment) by ::

  pip install fortnet-python

or into the user space issueing::

  pip install --user fortnet-python

Locally from Source
-------------------

Alternatively, you can install it locally from source, i.e. from the root folder
of the project::

  python -m pip install .

Testing
=======

The regression testsuite utilizes the `pytest` framework and may be executed by
::

  python -m pytest --basetemp=Testing
