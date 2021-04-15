.. _sec-fnetdata:
.. highlight:: none

####################
Generating a Dataset
####################

[Input: `recipes/fnetdata/`]

This chapter should serve as a tutorial guiding you through your first dataset
creation using the ``Fortformat`` Python class. A distinction is made between
two different types of target values; with separate sections being dedicated to:

* global system properties
* and atomic properties

After this tutorial, you will be able to create a Fortnet compatible dataset,
based on the output files of your simulation package of choice (e.g. VASP).


********************************************
Fortformat: Basic Fortnet Input Format Class
********************************************

This Python package provides a basic class which implements the associated
Fortnet input file format. The input features and targets, stored in lists
of Numpy arrays, conveniently get dumped to disk as ``fnetdata.xml`` files.


Installation
============

Please note, that this package has been tested for Python 3.X support. It
additionally needs Numerical Python (the Numpy module).

Since the ``Fortformat`` class, among others, expects the so-called `Atoms`
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


*****************
Global Properties
*****************

[Input: `recipes/fnetdata/globalTargets/`]

If a training on system-wide, global properties (e.g. the total energy) is
desired, this section is an ideal introduction to generating a suitable dataset.

As an application example, the :math:`E`-:math:`V` scan of a primitive silicon
unitcell in the diamond phase is used. The calculations were carried out by the
famous Vienna Ab initio Simulation Package (VASP)
:cite:`vasp1,vasp2,vasp3,vasp4` for next neighbor distances in the interval
:math:`[2.00,3.50]\,\mathrm{Å}` and a stepsize of :math:`0.05\,\mathrm{Å}`. The
raw data is stored at `recipes/fnetdata/globalTargets/vaspdata/` in the form of
structure information (POSCAR) and simulation output (OUTCAR). The following
Python script shows one possible way to get a dataset containing the total
energy of the respective system, based on this raw data.
 
.. code-block:: python

  #!/usr/bin/env python3

  '''
  Application example of the Fortformat class, based on a dataset
  that provides global system properties as target values to fit.
  '''

  import os
  import numpy as np
  from fortformat import Fortformat
  from ase.io.vasp import read_vasp, read_vasp_out

  def main():
      '''Main driver routine.'''

      nndists = np.arange(2.00, 3.50 + 0.05, 0.05)

      inpaths = [os.path.join(os.getcwd(), 'vaspdata', entry)
		 for entry in sorted(os.listdir('vaspdata'))]
      outpaths = [os.path.join(os.getcwd(), 'dataset', 'nndist_{:.3f}'
			       .format(nndist)) for nndist in nndists]

      strucs = []
      energies = np.empty((len(inpaths), 1))

      for ii, inpath in enumerate(inpaths):
	  os.makedirs(outpaths[ii], exist_ok=True)
	  struc = read_vasp(os.path.join(inpath, 'POSCAR'))
	  strucs.append(struc)
	  props = read_vasp_out(os.path.join(inpath, 'OUTCAR'))
	  energies[ii, 0] = props.get_total_energy()

      fnetdata = Fortformat(strucs, outpaths, targets=energies,
                            atomic=False, frac=True)
      fnetdata.dump()

  if __name__ == '__main__':
      main()

Following the necessary imports, the main method first generates the
corresponding next neighbor distances as already mentioned above. Two simple
list comprehensions further establish lists with the in- and output paths. While
iterating over all input paths, each corresponding output folder gets created
and the ASE `Atoms` object appended to an empty list of structures. The
individual total energies of the datapoints are stored in an empty Numpy array,
where the number of rows being determined by the number of datapoints and the
columns by the number of global targets per datapoint. Finally, a ``Fortformat``
object gets instantiated using the gathered informations, as well as providing
keyword arguments to determine if atomic properties are present (default: False)
and whether the coordinates should be saved in fractional or absolute format
(default: False).


*****************
Atomic Properties
*****************

[Input: `recipes/fnetdata/atomicTargets/`]

If training on atom specific properties (e.g. atomic forces) is desired, then
this section is an ideal introduction to generating a suitable dataset.

As an application example, the :math:`E`-:math:`V` scan of a primitive silicon
unitcell in the diamond phase is used. The calculations were carried out by the
famous Vienna Ab initio Simulation Package (VASP)
:cite:`vasp1,vasp2,vasp3,vasp4` for next neighbor distances in the interval
:math:`[2.00,3.50]\,\mathrm{Å}` and a stepsize of :math:`0.05\,\mathrm{Å}`. The
raw data is stored at `recipes/fnetdata/atomicTargets/vaspdata/` in the form of
structure information (POSCAR) and simulation output (OUTCAR). The following
Python script shows one possible way to get a dataset containing the total
energy per atom of the respective system, based on this raw data. Please note
that this is for demonstration purposes only and has no direct physical
relevance. A more sensible dataset could, for example, contain the atomic forces
as targets.

.. code-block:: python

  #!/usr/bin/env python3

  '''
  Application example of the Fortformat class, based on a dataset
  that provides atomic system properties as target values to fit.
  '''

  import os
  import numpy as np
  from fortformat import Fortformat
  from ase.io.vasp import read_vasp, read_vasp_out

  def main():
      '''Main driver routine.'''

      nndists = np.arange(2.00, 3.50 + 0.05, 0.05)

      inpaths = [os.path.join(os.getcwd(), 'vaspdata', entry)
		 for entry in sorted(os.listdir('vaspdata'))]
      outpaths = [os.path.join(os.getcwd(), 'dataset', 'nndist_{:.3f}'
			       .format(nndist)) for nndist in nndists]

      strucs = []
      energies = []

      for ii, inpath in enumerate(inpaths):
	  os.makedirs(outpaths[ii], exist_ok=True)
	  struc = read_vasp(os.path.join(inpath, 'POSCAR'))
	  strucs.append(struc)
	  props = read_vasp_out(os.path.join(inpath, 'OUTCAR'))
	  tmp = np.empty((len(struc), 1))
	  tmp[:, 0] = props.get_total_energy() / 2.0
	  energies.append(tmp)

      fnetdata = Fortformat(strucs, outpaths, targets=energies,
                            atomic=True, frac=True)
      fnetdata.dump()

  if __name__ == '__main__':
      main()

The procedure is nearly analogous to the global target example above: Following
the necessary imports, the main method first generates the corresponding next
neighbor distances as already mentioned above. Two simple list comprehensions
further establish lists with the in- and output paths. While iterating over all
input paths, each corresponding output folder gets created and the ASE `Atoms`
object appended to an empty list of structures. Since each of those structures
will in general have a different number of atoms, the target values are stored
in a list of Numpy arrays, where the number of rows being determined by the
number of atoms and the columns by the number of targets per atom. Finally, a
``Fortformat`` object gets instantiated using the gathered informations, as well
as providing keyword arguments to determine if atomic properties are present
(default: False) and whether the coordinates should be saved in fractional or
absolute format (default: False).
