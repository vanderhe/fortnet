.. _sec-fnetdata:
.. highlight:: none

##############################
Fnetdata: Generating a Dataset
##############################

[Input: `recipes/fortformat/fnetdata/`]

This chapter should serve as a tutorial guiding you through your first dataset
creation using the ``Fnetdata`` Python class. A distinction is made between
two different types of target values; with separate sections being dedicated to:

* global system properties
* and atomic properties
* a mixture of global and atomic properties

After this tutorial, you will be able to create a Fortnet compatible dataset,
based on the output files of your simulation package of choice (e.g. VASP).

.. _sec-fnetdata_globalTargets:

*****************
Global Properties
*****************

[Input: `recipes/fortformat/fnetdata/globalTargets/`]

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
  Application example of the Fnetdata class, based on a dataset
  that provides global system properties as target values to fit.
  '''

  import os
  import numpy as np
  from fortformat import Fnetdata
  from ase.io.vasp import read_vasp, read_vasp_out

  def main():
      '''Main driver routine.'''

      nndists = np.arange(2.00, 3.50 + 0.05, 0.05)

      inpaths = [os.path.join(os.getcwd(), 'vaspdata', entry)
		 for entry in sorted(os.listdir('vaspdata'))]

      strucs = []
      energies = np.empty((len(inpaths), 1))

      for ii, inpath in enumerate(inpaths):
	  strucs.append(read_vasp(os.path.join(inpath, 'POSCAR')))
	  props = read_vasp_out(os.path.join(inpath, 'OUTCAR'))
	  energies[ii, 0] = props.get_total_energy()

      fnetdata = Fnetdata(atoms=strucs, globaltargets=energies)
      fnetdata.dump('fnetdata.hdf5')

  if __name__ == '__main__':
      main()

Following the necessary imports, the main method first generates the
corresponding next neighbor distances as already mentioned above. A simple
list comprehension further establishes a list containing the input paths. While
iterating over all input paths, an ASE `Atoms` object gets appended to an empty
list of structures. The individual total energies of the datapoints are stored
in an empty Numpy array, where the number of rows being determined by the number
of datapoints and the columns by the number of global targets per datapoint.
Finally, an ``Fnetdata`` object gets instantiated using the gathered
information.


*****************
Atomic Properties
*****************

[Input: `recipes/fortformat/fnetdata/atomicTargets/`]

If training on atom specific properties (e.g. atomic forces or charges) is
desired, then this section is an ideal introduction to generating a suitable
dataset.

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
or charges as targets.

.. code-block:: python

  #!/usr/bin/env python3

  '''
  Application example of the Fnetdata class, based on a dataset
  that provides atomic system properties as target values to fit.
  '''

  import os
  import numpy as np
  from fortformat import Fnetdata
  from ase.io.vasp import read_vasp, read_vasp_out

  def main():
      '''Main driver routine.'''

      nndists = np.arange(2.00, 3.50 + 0.05, 0.05)

      inpaths = [os.path.join(os.getcwd(), 'vaspdata', entry)
		 for entry in sorted(os.listdir('vaspdata'))]

      strucs = []
      energies = []

      for ii, inpath in enumerate(inpaths):
	  struc = read_vasp(os.path.join(inpath, 'POSCAR'))
	  strucs.append(struc)
	  props = read_vasp_out(os.path.join(inpath, 'OUTCAR'))
	  tmp = np.empty((len(struc), 1))
	  tmp[:, 0] = props.get_total_energy() / 2.0
	  energies.append(tmp)

      fnetdata = Fnetdata(atoms=strucs, atomictargets=energies)
      fnetdata.dump('fnetdata.hdf5')

  if __name__ == '__main__':
      main()

The procedure is nearly analogous to the global target example above: Following
the necessary imports, the main method first generates the corresponding next
neighbor distances as already mentioned above. A simple list comprehension
further establishes a list containing the input paths. While iterating over all
input paths, an ASE `Atoms` object gets appended to an empty list of structures.
Since each of those structures will in general have a different number of atoms,
the target values are stored in a list of Numpy arrays, where the number of rows
being determined by the number of atoms and the columns by the number of targets
per atom. Finally, an ``Fnetdata`` object gets instantiated using the gathered
information.


********************
Weighting Datapoints
********************

[Input: `recipes/fortformat/fnetdata/weighting/datapoints/`]

There are conceivable situations in which weighting individual datapoints makes
sense. The detour via the increased insertion of a datapoint is not only
cumbersome but also inefficient, since exactly the same input features
(e.g. ACSF) and gradients would be calculated multiple times. To elegantly
circumvent this, ``Fnetdata`` and ``Fortnet`` offer the possibility of
individually weighting certain datapoints of a dataset. After a Fortformat
object has been instantiated, the desired weights can be handed over via a
setter function. The following code snippet shows what this could look like:

.. code-block:: python

  # start with homogeneous weighting
  weights = np.ones((31,), dtype=int)
  # possibly, certain datapoints are more important
  weights[4:13] = 3

  fnetdata = Fnetdata(atoms=strucs, globaltargets=energies)
  fnetdata.weights = weights
  fnetdata.dump('fnetdata.hdf5')

For Fortformat to correctly recognize the weights, they must be specified as a
onedimensional list or Numpy array of positive integers. If these requirements
are not met, an error message is issued, so that nothing can terribly go wrong
(fingers crossed).


**************************
Weighting Atomic Gradients
**************************

[Input: `recipes/fortformat/fnetdata/weighting/gradients/`]

Further, there might be a need for different weighting of atomic contributions
in the training process. This allows to change the contribution of specific
atoms to the training process, as well as to completely `switch off` atoms, if
the respective target would not be defined. Therefore, ``Fnetdata`` and
``Fortnet`` offer the possibility of setting atom-resolved weights after a
Fortformat object has been instantiated. The desired weights can be handed over
via a setter function. The following code snippet shows what this could look
like:

.. code-block:: python

  # fix random seed for reproduction purposes
  np.random.seed(42)

  atomicweights = []
    .
    .
  for ii, atom in enumerate(atoms):
          .
	  .
      # float-valued atomic gradient weighting in interval [1, 10]
      atomicweights.append(np.asfarray(
	  np.random.randint(1, 10, len(atom), dtype=int)))

  fnetdata = Fnetdata(atoms=strucs, globaltargets=energies)
  fnetdata.atomicweights = atomicweights
  fnetdata.dump('fnetdata.hdf5')

Alternatively, the weights can be boolean-valued. This allows the contributions
of individual atoms to be switched `on` or `off`. Currently, these are
internally converted to floats 0.0 (False) and 1.0 (True), so there is no
performance advantage. In the future, however, the corresponding gradient
calculations will be skipped and thus a significant performance gain achieved:

.. code-block:: python

  # fix random seed for reproduction purposes
  np.random.seed(42)
  sample = [True, False]

  atomicweights = []
    .
    .
  for ii, atom in enumerate(atoms):
          .
	  .
      # randomly activate/deactivate atomic contributions
      atomicweights.append(np.random.choice(sample, size=len(atom)))

  fnetdata = Fnetdata(atoms=strucs, globaltargets=energies)
  fnetdata.atomicweights = atomicweights
  fnetdata.dump('fnetdata.hdf5')

For Fortformat to correctly recognize the weights, they must be specified as a
onedimensional list of lists/numpy arrays of values :math:`\geq` 0. If these
requirements are not met, an error message is issued, so that nothing can
terribly go wrong (again, fingers crossed).

.. _sec-fnetdata_extFeatures:

************************
External Atomic Features
************************

[Input: `recipes/fortformat/fnetdata/extfeatures/`]

Since currently only the Atom-Centered Symmetry Functions are available as a
mapping of the geometries to infer suitable network inputs, ``Fnetdata``, as
well as ``Fortnet``, offer the possibility of processing user-specified external
atomic features. Thus every kind of imaginable input features can be used in the
training and prediction process, which significantly expands the versatility of
Fortnet. Of course, the user is responsible for checking the suitability of the
features handed over. Charge specifications like the Mulliken populations of the
individual atoms are conceivable.

The transfer of the selected features to the Fortformat class is
straightforward. The keyword argument ``features`` expects a list of Numpy
arrays, where the first dimension corresponds to the number of atoms of the
associated geometry and the second to the number of features per atom. The
example below shows how such a specification could look like. Random numbers are
used as features, which therefore only serve a demonstrative purpose.

.. code-block:: python

  #!/usr/bin/env python3

  '''
  Application example of the Fnetdata class, based on a dataset
  that provides global system properties as target values to fit.
  The dataset is extended by user specified external atomic features.
  '''

  import os
  import numpy as np
  from fortformat import Fnetdata
  from ase.io.vasp import read_vasp, read_vasp_out

  def main():
      '''Main driver routine.'''

      np.random.seed(42)

      inpaths = [os.path.join(os.getcwd(), '../globalTargets/vaspdata', entry)
		 for entry in sorted(os.listdir('../globalTargets/vaspdata'))]

      strucs = []
      features = []
      energies = np.empty((len(inpaths), 1))

      for ii, inpath in enumerate(inpaths):
	  struc = read_vasp(os.path.join(inpath, 'POSCAR'))
	  strucs.append(struc)
	  props = read_vasp_out(os.path.join(inpath, 'OUTCAR'))
	  energies[ii, 0] = props.get_total_energy()
	  features.append(np.random.random_sample((len(struc), 3)))

      fnetdata = Fnetdata(atoms=strucs, globaltargets=energies,
                          features=features)
      fnetdata.dump('fnetdata.hdf5')

  if __name__ == '__main__':
      main()
