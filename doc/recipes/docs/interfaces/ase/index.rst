.. _sec-interfaces-ase:

###################################
Atomic Simulation Environment - ASE
###################################

The Atomic Simulation Environment - ASE is a set of Python based tools and 
modules for setting up, manipulating, running, visualizing and analyzing 
atomistic simulations (cf. `ASE documentation 
<https://wiki.fysik.dtu.dk/ase/>`_). 

Currently, file-IO based communication between ASE and Fortnet is available.
Further information can be found in the section linked below.

.. toctree::
   :maxdepth: 1

   fileio/fileio.rst

Note: Before going through the following sections, please make sure that you
have installed a working version of the `ASE` and `fortnet-python` package. If
you are wondering how to
`install ASE <https://wiki.fysik.dtu.dk/ase/install.html>`_
or `install fortnet-python <https://github.com/vanderhe/fortnet-python/blob/
master/README.rst>`_, please consult the corresponding documentation.

Installation
============

Please note, that the corresponding package has been tested for Python 3.X
support. Its usage additionally requires

  - `Atomic Simulation Environment
    <https://wiki.fysik.dtu.dk/ase/install.html>`_ (`ase`)
  - `Fortnet Python tools <https://github.com/vanderhe/fortnet-python>`_
    (`fortnet-python`)

as well as the `pytest <https://docs.pytest.org/en/6.2.x/>`_ framework in order
to run the regression tests.

Via the Python Package Index
----------------------------

First, make sure you have an up-to-date version of pip installed::

  python -m pip install --upgrade pip

The package can be downloaded and installed via pip into the active Python
interpreter (preferably using a virtual python environment) by ::

  pip install fortnet-ase

or into the user space issueing::

  pip install --user fortnet-ase

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
