.. _sec-firstpredict:
.. highlight:: none

******************************
First Predictions with Fortnet
******************************

[Input: `recipes/basics/firstpredict/`]

This chapter should serve as a tutorial guiding you through your first
predictions using Fortnet. The network from the
:ref:`previous section <sec-firsttrain>`, trained on the :math:`E`-:math:`V`
scan of a primitive silicon unitcell in the diamond phase, is used as a starting
point. The procedure is split into three major steps we are already familiar
with:

* telling Fortnet what to do,
* running Fortnet,
* analysing the results.

Providing the Input
===================

[Input: `recipes/basics/firstpredict/validate/`]

Fortnet accepts the input in the Human-readable Structured Data (HSD) format.
The input file must be called `fortnet_in.hsd`.  The input file used in this
example looks as follows::

  Data {
    Dataset = fnetdata.hdf5
    NetstatFile = fortnet.hdf5
  }

  Options {
    ReadNetStats = Yes
    Mode = validate
  }

The order of the specified blocks in the HSD input is arbitrary. You are free to
capitalise the keywords as you like, since they are case-insensitive. This is
not valid however for string values, especially if they are specifying file
names.

So let's have a look at the two necessary input blocks, ``Data`` and
``Options``.

Data
----
::

  Data {
    Dataset = fnetdata.hdf5
    NetstatFile = fortnet.hdf5
  }

The ``Data`` block looks exactly the same as for the training process, because
we first want to check that the network has actually learned something and is
able to reproduce the training dataset with sufficient accuracy. With other
words, we are validating the resulting network potential, which leads us to the
running mode specified in the next HSD block below.

Options
-------
::

  Options {
    ReadNetStats = Yes
    Mode = validate
  }

The basic program behavior gets defined in the ``Option`` block of the input,
starting with the running mode of Fortnet. Since we want to predict structures
based on an existing network potential, we have to set the ``ReadNetStats``
entry to `Yes` in order to read in the netstat file. Subsequently the ``Mode``
is set to `validate`, so Fortnet will do a prediction run and additionally list
the corresponding target values in the ``fnetout.hdf5`` output file for later
comparison.

.. note::

   When ``ReadNetStats`` is set to `No` in validation or prediction mode,
   Fortnet will issue a warning that tells the user that this isn't a valid
   combination and overwrites the input to `Yes`.


Running Fortnet
===============
At this point you are ready to execute Fortnet. To do so, invoke the ``fnet``
binary without any arguments in the directory containing the ``fortnet_in.hsd``
file. As mentioned above, Fortnet writes some information to the standard
output. Therefore it is recommended to tee this output for later investigation::

  fnet | tee output

In most cases Fornet will be compiled with MPI parallelism enabled. To make use
of the associated speedup, issue::

  mpirun -np 4 fnet | tee output

or something equivalent. Note: It may be necessary to provide the absolute path
to the ``fnet`` binary in this case.


Examining the Output
====================
Fortnet uses two output channels: 1) the standard output (which you should
redirect into a file to keep for later evaluation) and 2) various output files.
Below, these two channels are routhly explained for the prediction scenario.

Standard Output
---------------
Most of the standard output is identical to what was described in the previous
section. In the ``Initialisation`` section, however, the entry
`read initial netstats` will be set to `True` or `T` respectively.
::

  Initialisation

  running in validation mode
  random seed: 571070859
  read initial netstats: T

Also, instead of the output of the training process, you will now see the simple
message `Start feeding...` that indicates the start of the feeding process.
::

  Start feeding...done

When finished, a `done` will be appended and the predictions written to disk
(c.f. next section).

Fnetout
-------
The ``fnetout.hdf5`` file is the most important output of Fortnet as it contains
all the predictions made. In validation mode this file will also contain the
target values provided by the dataset, whereas in prediction mode thoose exact
values are generally unknown and therefore not contained in the output. Again,
feel free to open the HDF5 file with your viewer of choice. The following script
shows how to extract the predictions and targets from the output file by using
the ``Fortformat`` Python package that ships with Fortnet:

.. code-block:: python

  #!/usr/bin/env python3

  '''
  Application example of the Fortformat package, based on an output
  file that contains network predictions and corresponding targets.
  '''

  import numpy as np
  import matplotlib.pyplot as plt
  from fortformat import Fnetout

  def main():
      nndists = np.arange(2.10, 3.30 + 0.05, 0.05)

      fnetout = Fnetout('fnetout.hdf5')
      globalpredictions = fnetout.globalpredictions
      globaltargets = fnetout.globaltargets

      plt.figure(figsize=(7, 5))
      plt.title('Comparison of Neural Network Predictions with Targets')
      plt.xlabel(r'Nearest Neighbour Distance [$\mathrm{\AA}$]')
      plt.ylabel('Total Energy [eV / Atom]')

      plt.plot(nndists, globalpredictions / 2.0, color='blue', label='NN')
      plt.scatter(nndists, globaltargets / 2.0, s=10, color='black', label='DFT')

      plt.tight_layout()
      plt.legend()
      plt.savefig('comparison.svg', dpi=900, format='svg')

  if __name__ == '__main__':
      main()

If the predictions and targets are being
plotted, an excellent agreement will be observed:

.. figure:: ../_figures/basics/firstpredict/energy-volume-scan.svg
   :width: 100%
   :align: center
   :alt: Comparison of neural network predictions with targets.

As a further analysis, the energies of next neighbor distances beyond the
training interval can be predicted. To do so, we finally got to use the pure
prediction mode by setting ``Mode`` of the ``Option`` block to `predict`. The
corresponding figure below impressively shows a major weakness of neural
networks, their poor extrapolation capabilities:

[Input: `recipes/basics/firstpredict/predict/`]

.. figure:: ../_figures/basics/firstpredict/e-v-scan_plus_extrap.svg
   :width: 100%
   :align: center
   :alt: Comparison of neural network predictions with targets.

Outside the next neighbor distances for which there was available data in the
training process (visualized by the vertical, dashed lines), there is a
significant deviation between the predictions and reference values. This is
something that must always be considered when dealing with neural networks.
