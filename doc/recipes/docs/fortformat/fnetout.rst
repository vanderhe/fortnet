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
datapoints the network was trained on, the type of training targets (atomic or
global), the predictions of the network potential as well as corresponding
targets if provided (only for validation mode).

The following Python script shows how to extract the aforementioned information
based on the :math:`E`-:math:`V` scan example of a primitive silicon unitcell:

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
      print('Running mode: ', mode)

      ndatapoints = fnetout.ndatapoints
      print('Number of datapoints in training: ', ndatapoints)

      targettype = fnetout.targettype
      print('Type of targets: ', targettype)

      predictions = fnetout.predictions
      print("Fortnet's predictions: ", predictions)

      targets = fnetout.targets
      print('Targets while trained: ', targets)

  if __name__ == '__main__':
      main()
