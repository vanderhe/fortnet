.. _sec-firsttrain:
.. highlight:: none

***************************
First Training with Fortnet
***************************

[Input: `recipes/basics/firsttrain/`]

This chapter should serve as a tutorial guiding you through your first network
optimization using Fortnet. As an exemplary dataset, the :math:`E`-:math:`V`
scan of a primitive silicon unitcell in the diamond phase is used. The procedure
is split into three major steps:

* providing an appropriate input to Fortnet,
* actually running Fortnet,
* finally, analysing the results

After this tutorial, you will therefore have already become familiar with all of
the basic features of Fortnet and subsequently start your own project.

Providing the Input
===================

Fortnet accepts the input in the Human-readable Structured Data (HSD) format.
The input file must be called `fortnet_in.hsd` and in this example looks as
follows::

  Network = BPNN {
    Hidden = 2 2
    Activation = tanh
  }

  Features {
    Mapping = ACSF {
      Reduce = Yes
      Standardization = Yes
      Function = Auto {
	RCut = 4.0
	NRadial = 5
	NAngular = 4
      }
    }
  }

  Training = LBFGS {
    Threshold = 1e-08
    NIterations = 5000
    NPrintout = 1000
    NSaveNet = 1000
    MinDisplacement = 1e-10
    MaxDisplacement = 5e-02
    LineMin = Yes
    Memory = 1000
    Loss = mse
  }

  Data {
    Dataset = fnetdata.hdf5
    NetstatFile = 'fortnet.hdf5'
  }

  Options {
    Mode = train
    ReadNetStats = No
    RandomSeed = 123456
  }

The order of the specified blocks in the HSD input is arbitrary. You are free to
capitalise the keywords as you like, since they are case-insensitive. This is
not valid however for string values, especially if they are specifying file
names. Furthermore, it is possible to put arbitrary comments in the HSD input
after a hash-mark (``#``) character. Everything between this character and the
end of the current line is ignored by the parser.

So let's have a look at the input blocks, one by one.

Network
-------
::

  Network = BPNN {
    Hidden = 2 2
    Activation = tanh
  }

The ``Network`` block specifies the neural network architecture to use.
Currently, only the Behler-Parrinello-Neural-Network (BPNN) :cite:`bpnn` type
is implemented. It is assumed that all sub-nn's have the same internal
structure, i.e. the same number of hidden layers and neurons per layer. The
``Hidden`` child node controls said parameters by expecting a list of positive
integer values, where each value corresponds to a hidden layer with the
specified number of neurons. To determine a neuron status, the activation or
transfer function is essential. Its type is controlled by the ``Activation``
entry. In this case, let's use the hyperbolic tangent. For complete list of
activation functions, please consult the corresponding
:ref:`section <sec-transfer>`.

Features
--------

::

  Features {
    Mapping = ACSF {
      Reduce = Yes
      Standardization = Yes
      Function = Auto {
	RCut = 4.0
	NRadial = 5
	NAngular = 4
      }
    }
  }

Fortnet tries to infer physical or chemical properties of your systems based on
structural information, i.e. the atom types and coordinates. Since these raw
values are unsuitable as network inputs, for several reasons, they have to get
mapped to translational, rotational and commutation (same type) invariant
values. One famous set of functions that fulfills this purpose are the so-called
Atom-centered symmetry functions (ACSF) by J. Behler :cite:`acsf`. Fortnet
currently implements radial :math:`G_1, G_2, G_3` and angular :math:`G_4, G_5`
functions, as denoted in the original ACSF paper. In this case Fortnet's
automatic parameter generation scheme is used to achieve a decent coverage of
the cutoff sphere by utilizing :math:`G_2` and :math:`G_5` functions. Therefore,
only the number of radial (``NRadial``) and angular (``NAngular``), as well as
the cutoff radius (``RCut``), needs to be specified. The unit of the cutoff
radius is Angstrom. Due to the nature of the ACSF it is likely to get input
values of very different orders of magnitude. To compensate for this and achieve
an improvement in convergency and overall stability, it is possible to apply a
simple z-score standardization in the background, before feeding the network.
This behavior is controlled via the ``Standardization`` option. The ``Reduce``
entry determines whether the ACSF functions should be element resolved or
unresolved. In the latter case (``Reduce`` = Yes) the calculated neighbor lists
would contain all the atoms regardless of their type which leads to a
significant reduction of the input features. However, this would require the
weighting of individual summands with :ref:`atomic prefactors <sec-acsf_atomid>`
since otherwise contradictory input features would arise. Since the dataset at
hand only contains silicon atoms, this parameter may be ignored for now.


Training
--------
::

  Training = LBFGS {
    Threshold = 1e-08
    NIterations = 5000
    NPrintout = 1000
    NSaveNet = 1000
    MinDisplacement = 1e-10
    MaxDisplacement = 5e-02
    LineMin = Yes
    Memory = 1000
    Loss = mse
  }

To successively optimize the weight and bias network parameters during the
training iterations, Fortnet provides different algorithms. In this example
a limited memory implementation of the Broyden–Fletcher–Goldfarb–Shanno
algorithm (L-BFGS) is used. For a complete list of the available optimizers,
please consult the corresponding :ref:`optimizer <sec-optimizer>` section. Every
optimizer provides two options to controll when to end the training process, the
``Threshold`` and maximum number of iterations (``NIterations``). The training
will be terminated as soon as one of the conditions is fulfilled. Furthermore,
the number of training iterations must be specified, after which the current
loss value and gradient gets printed to stdout (``NPrintout``) and the current
network status is written out (``NSaveNet``). For a list of available loss
functions, consult the dedicated :ref:`Loss Functions <sec-loss>` section. The
remaining settings of the example above are optional and described in the
corresponding :ref:`L-BFGS <sec-optimizer>` optimizer subsection.


Data
----
::

  Data {
    Dataset = fnetdata.hdf5
    NetstatFile = fortnet.hdf5
  }

Since the provision of high quality data is key when dealing with neural
networks in general, let's have a look at the data block and how to hand over a
dataset. Most important, the ``Dataset`` entry must be a string pointing to a
compatible HDF5 dataset file (in this case ``fnetdata.hdf5``). A fundamental
design decision of Fortnet is not to provide native support for the output files
of popular simulation packages directly. Instead, a separate input format is
used and a corresponding Python class is provided which, based on the
Atomic Simulation Environment (`ASE <https://wiki.fysik.dtu.dk/ase/>`_) that is
also implemented in Python, enables a dataset to be generated easily. To see how
you get from the output files of your simulation package of choice to a Fortnet
compatible dataset, please consult the
:ref:`Fnetdata: Generating a Dataset <sec-fnetdata>` section.

Another useful feature is that the loss function of an external validation
dataset, that is not included in the optimization prozess, can be monitored
during training. To utilize this so-called validation-monitoring, e.g. for early
stopping purposes, provide an additional pathfile via the ``Validset`` entry::

  Data {
       .
    Validset = fnetvdata.hdf5
  }

In this case a dataset file named `fnetdata.hdf5` is present in the same folder
as the ``fortnet_in.hsd`` input. Feel free to have a look at its content by
using your HDF5 viewer of choice.

In addition, the ``Data`` block also handles the filename of the so-called
`netstat` files of the Fortnet world. They define the whole network status and
will be needed for a later restart of the training process or predictions based
on the created potential.


Options
-------
::

  Options {
    Mode = train
    ReadNetStats = No
    RandomSeed = 123456
  }

The basic program behavior gets defined in the ``Option`` block of the input,
starting with the running mode of Fortnet. There are three valid options:
`train`, `validate`, `predict`. As in this example, the `train` mode will
optimize the network with respect to the targets provided by the dataset. A
resumption of the training process based on existing `netstat` file would be
requested by setting the ``ReadNetStats`` entry to `Yes`. To validate the
resulting networks or to predict structures with unknown properties, the
other two modes are used and explained in the
:ref:`First Predictions with Fortnet <sec-firstpredict>` section.

The reproducibility of results is particularly important in scientific fields of
application. To meet this requirement, Fortnet provides a ``RandomSeed`` entry.
By setting a seed you define the initial state of the luxury random number
generator :cite:`ranlux1,ranlux2,ranlux3` that is working in the background and
is responsible for the outcome of the initialization of the sub-nn's and
therefore the training process in general. This is an optional entry and
randomly generated if not set by the user. Since Fortnet prints out the random
seed of the current run you may need this for later reproduction of results.

.. warning::
   A few warning words about the reproducibility: In theory all the results you
   obtain using Fortnet are reproducible since the ``RandomSeed`` entry enables
   the user to define the initial state of the random number generators used by
   the project. However, due to the non-commutativity of floating-point
   operations it has been observed that reproducibility is given for a fixed
   machine, compiler and number of MPI-processes, but as soon as one of these
   parameters changes you will get different results.


Running Fortnet
===============
As soon as all files have been generated and are present in their correct
location, you are ready to execute Fortnet. To do so, invoke the ``fnet`` binary
without any arguments in the directory containing the ``fortnet_in.hsd`` file.
As mentioned above, Fortnet writes some information to the standard output.
Therefore it is recommended to tee this output for later investigation::

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
These two channels will now be outlined, within the context of a training
scenario, below.

Standard Output
---------------
In the following, the standard output, gets broken down and explained piece by
piece, in the order as it appears on the screen, starting with the header::

  |==============================================================================|
  |  Fortnet - A BPNN Implementation, Version 0.3                                |
  |                                                                              |
  |  Copyright (C) 2020 - 2021  T. W. van der Heide                              |
  |==============================================================================|

  date: 15.08.2021
  time: 13:13:08, +0200

As you may have seen, nothing spectacular is happening here. Nevertheless, the
version number as well as date and time of the binary execution can be important
information in retrospect.

::

  Interpreting input file 'fortnet_in.hsd'
  Checking Input Consistency...passed
  Processed input written as HSD to 'fortnet_pin.hsd'

  --------------------------------------------------------------------------------

As the next step, Fortnet parses and interprets the ``fortnet_in.hsd`` input
file and carries out some basic consistency checks on the obtained parameters.
Additionally the input as Fortnet sees and interprets it gets stored in the
``fortnet_pin.hsd`` file.

You will also see a list of information from the HSD input, as printed below::

  Initialisation

  running in training mode
  random seed: 123456
  read initial netstats: F

  --------------------------------------------------------------------------------

  Sub-NN Details

  inputs: 9
  hidden layers: 2 2
  outputs: 1

  activation: tanh

  --------------------------------------------------------------------------------

  ACSF Mappings

  species-resolved: F

  nr. of radial functions: 5
  nr. of angular functions: 4

  g2: rc = 7.558904, rs = .000000, eta = .805987,
      atomId = 0
  g2: rc = 7.558904, rs = 1.889726, eta = .805987,
      atomId = 0
  g2: rc = 7.558904, rs = 3.779452, eta = .805987,
      atomId = 0
  g2: rc = 7.558904, rs = 5.669178, eta = .805987,
      atomId = 0
  g2: rc = 7.558904, rs = 7.558904, eta = .805987,
      atomId = 0
  g5: rc = 7.558904, lambda = 1.000000, eta = .080599, xi = 1.000000
      atomId = 0
  g5: rc = 7.558904, lambda = -1.000000, eta = .080599, xi = 1.000000
      atomId = 0
  g5: rc = 7.558904, lambda = 1.000000, eta = .080599, xi = 16.000000
      atomId = 0
  g5: rc = 7.558904, lambda = -1.000000, eta = .080599, xi = 16.000000
      atomId = 0
  --------------------------------------------------------------------------------

  Dataset Information

  found: 25 datapoints (25 unique ones)
  in file: fnetdata.hdf5
  total sub-nn parameters: 29
  targets per parameter: .8621

  --------------------------------------------------------------------------------

The entry ``targets per parameter`` is of particular importance. Based on this
ratio you can roughly deduce whether the selected network size is suitable
regarding the dataset that was provided. It is calculated in terms of unique
datapoints, by solely considering the unweighted geometry-target pairs.

Up to this stage of binary execution, the input was parsed and the dataset read.
The ``Calculating ACSF`` statement tells us, that Fortnet has started to map the
structure information to input-suitable ACSF values. As soon as the word `done`
appears, this process is complete and the training process starts::

  Calculating ACSF...done
  Starting training...

       iTrain           MSE-Loss          Gradients
  --------------------------------------------------------------------
	1000        0.187303E-04       0.548379E-03
	2000        0.224850E-04       0.121142E-02
	3000        0.411225E-05       0.315495E-03
	4000        0.118531E-05       0.820748E-03
	5000        0.142334E-06       0.205431E-03
  --------------------------------------------------------------------

  Training finished (max. Iterations reached)

  --------------------------------------------------------------------

  Loss Analysis (global min.)

  iTrain: 5000, Loss: 1.423336E-07

  --------------------------------------------------------------------

While the training process is running, the trajectory of the loss function and
the total gradient of the network parameters are printed regularly, depending on
the ``NPrintout`` setting of the ``Training`` block. In this case, the
termination criterion is the maximum number of training iterations. After
completion of the training, the iteration with the lowest loss value is written
out. 

Output Files
------------
Depending on the setting of the program behavior in the input file (i.e. running
mode), different output files are created. Running the current example there
will be a single file written to disk, appart from the redirected standard
output: `fortnet.hdf5`. The average user does not have to look into this file.
It solely contains information regarding the program state, which are necessary
for a later resumption of the training process or for predictions based on the
resulting network potential.

In fact, the relevant output ``fnetout.hdf5`` is only created in validation or
prediction mode and introduced in the :ref:`next section <sec-firstpredict>`.

If the total trajectory of the loss function and total gradient is of interest,
it can be written out as ``iterout.dat`` by setting the corresponding entry
(default: No)::

  Options {
       .
       .
       .
    WriteIterationTrajectory = Yes
  }

The column order of the output in ``iterout.dat`` is analogous to the standard
output.
