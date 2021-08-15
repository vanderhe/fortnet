*******************************************
FORTFORMAT: Basic Fortnet IO Format Classes
*******************************************

These basic Python classes implement the Fortnet input and output file format.
The Fnetdata class enables to create compatible HDF5 datasets, whereas the
Fnetout class extracts certain properties of the HDF5 output for later analysis.
The corresponding Python classes are referred to as:

Fnetdata
  Basic Fortnet input format Python class.

Fnetout
  Basic Fortnet output format Python class.


Installation
============

Please note, that this package has been tested for **Python 3.X**
support. It additionally needs Numerical Python (the Numpy module).

System install
--------------

You can install the script package via the standard 'python setup'
mechanism. If you want to install it system-wide into your normal
python installation, you simply issue::

  python setup.py install

with an appropriate level of permission.

Local install
-------------

Alternatively, you can install it locally in your home space, e.g.::

  python setup.py install --user

If the local python install directory is not in your path, you should
add this. For the bash shell you should include in .bashrc::

  export PATH=$PATH:/home/user/.local/bin


For developers
--------------

To perform pylint static checking from the top level directory of the
Fortnet project, use

pylint3 --rcfile utils/srccheck/pylint/pylintrc-3.ini tools/fortformat/src/*
