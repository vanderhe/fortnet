********************************************
FORTFORMAT: Basic Fortnet Input Format Class
********************************************

This Python package provides a basic class which implements the associated
Fortnet input file format. The input features and targets, stored in lists
of Numpy arrays, conveniently get dumped to disk as simple text files. The
corresponding Python class is named:

Fortformat
  Basic Fortnet input format Python class.


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
