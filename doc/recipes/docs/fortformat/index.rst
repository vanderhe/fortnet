.. _sec-fortformat:
.. highlight:: none

##########
Fortformat
##########

[Input: `recipes/fortformat/`]

This chapter should serve as a tutorial guiding you through your first dataset
creation and data extraction after running Fortnet by using two Python classes,
provided by the ``Fortformat`` Python package:

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

Please note, that this package has been tested for Python 3.X support. It
additionally needs numerical Python (the `Numpy` module) as well as `h5py`.

Since the ``Fnetdata`` class, among others, expects the so-called `Atoms`
objects of the Atomic Simulation Environment
(`ASE <https://wiki.fysik.dtu.dk/ase/>`_) as an input, sooner or later this
dependency will also have to be satisfied.


System install
--------------

You can install the script package via the standard Python setup mechanism. If
you want to install it system-wide into your normal Python installation, you
simply issue
::

  python setup.py install

at `tools/fortformat/` with an appropriate level of permission.

Local install
-------------

Alternatively, you can install it locally in your home space, e.g.::

  python setup.py install --user

If the local Python install directory is not in your path, you should add this.
For the bash shell you should include the following line in the .bashrc::

  export PATH=$PATH:/home/user/.local/bin
