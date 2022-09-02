.. _sec-fnetout:
.. highlight:: none

##############################
Fnetout: Extracting the Output
##############################

[Input: `recipes/fortformat/fnetout/`]

This chapter should serve as a tutorial guiding you through the process of
extracting results from the main output file ``fnetout.hdf5`` of Fortnet by
using the ``Fnetout`` Python class.

*********************
Extracting Properties
*********************

[Input: `recipes/fortformat/fnetout/`]

To fetch information from an ``fnetout.hdf5`` output file, the ``Fnetout`` class
provides several properties that may be extracted, including the mode of the
Fortnet run that produced the output file (predict or validate), the number of
datapoints the network was trained on, the number of system-wide targets (e.g.
total energies) the network was trained on, the number of atomic targets (e.g.
atomic forces) the network was trained on, atomic force predictions (if
applicable) and (atom-resolved) system-wide and atomic predictions of the
network potential, as well as corresponding targets if provided (only for
validation mode).

The following Python script shows how to extract the aforementioned information:

.. code-block:: python

  #!/usr/bin/env python3

  '''
  Application example of the Fnetout class, based on an output
  file that provides the network predictions and targets, relevant
  information are extracted and printed to the standard output.
  '''

  from fortformat import Fnetout

  def main():
      '''Main driver routine.'''

      fnetout = Fnetout('fnetout.hdf5')

      mode = fnetout.mode
      print('Running mode:\n', mode)

      ndatapoints = fnetout.ndatapoints
      print('Number of datapoints in dataset:\n', ndatapoints)

      nglobaltargets = fnetout.nglobaltargets
      print('Number of system-wide targets (e.g. total energies):\n',
	    nglobaltargets)

      natomictargets = fnetout.natomictargets
      print('Number of atomic targets (e.g. atomic forces):\n',
	    natomictargets)

      globaltargets = fnetout.globaltargets
      print('System-wide targets:\n', globaltargets)

      atomictargets = fnetout.atomictargets
      print('Atomic targets: ', atomictargets)

      tforces = fnetout.tforces
      print('Whether atomic forces are present:\n', tforces)

      forces = fnetout.forces
      print('Atomic forces:\n', forces)

      atomicpredictions = fnetout.atomicpredictions
      print("Fortnet's predictions of atomic targets:\n", atomicpredictions)

      globalpredictions = fnetout.globalpredictions
      print("Fortnet's predictions of system-wide targets:\n", globalpredictions)

      globalpredictions_atomic = fnetout.globalpredictions_atomic
      print("Fortnet's atom-resolved predictions of system-wide targets:\n",
	    globalpredictions_atomic)

  if __name__ == '__main__':
      main()
