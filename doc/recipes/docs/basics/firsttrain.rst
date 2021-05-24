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
    Activation = 'tanh'
  }

  Mapping = ACSF {
    NRadial = 5
    NAngular = 4
    RCut = 4.0
    Standardization = Yes
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
    Loss = 'rms'
  }

  Data {
    Dataset = 'training_data'
    Standardization = No
    NetstatFiles = Type2FileNames {
      Prefix = "./"
      Suffix = ".net"
      LowerCaseTypeName = No
    }
  }

  Options {
    Mode = 'train'
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
    Activation = 'tanh'
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

Mapping
-------

::

  Mapping = ACSF {
    NRadial = 5
    NAngular = 4
    RCut = 4.0
    Standardization = Yes
  }

Fortnet tries to infer physical or chemical properties of your systems based on
structural information, i.e. the atom types and coordinates. Since these raw
values are unsuitable as network inputs, for several reasons, they have to get
mapped to translational, rotational and commutation (same type) invariant
values. One famous set of functions that fulfills this purpose are the so-called
Atom-centered symmetry functions (ACSF) by J. Behler :cite:`acsf`. Fortnet
currently implements radial :math:`G_2` and angular :math:`G_5` functions, as
denoted in the original ACSF paper. Their respective parameters are calculated
automatically by Fortnet, so that a decent coverage of the sphere defined by
the cutoff radius is guaranteed. Therefore, only the number of radial
(``NRadial``) and angular (``NAngular``), as well as the cutoff radius
(``RCut``), needs to be specified. The unit of the cutoff radius is Angstrom.
Due to the nature of the ACSF it is likely to get input values of very different
magnitudes of order. To compensate for this and achieve an improvement in
convergency and overall stability, it is possible to apply a simple z-score
standardization in the background, before feeding the network. This behavior is
controlled via the ``Standardization`` option.


Training
--------
::

  Training = LBFGS {
    Threshold = 1e-08
    NIterations = 10000
    NPrintout = 10
    NSaveNet = 100
    MinDisplacement = 1e-10
    MaxDisplacement = 5e-02
    LineMin = Yes
    Memory = 1000
    Loss = 'rms'
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
    Dataset = 'training_data'
    Standardization = No
    NetstatFiles = Type2FileNames {
      Prefix = "./"
      Suffix = ".net"
      LowerCaseTypeName = No
    }
  }

Since the provision of high quality data is key when dealing with neural
networks in general, let's have a look at the data block and how to hand over a
dataset. Most important, the ``Dataset`` entry must be a string pointing to a
file that contains all the paths to the so called ``fnetdata.xml`` files. Each
of those files defines a datapoint that consists of a geometry and target
values to optimize the network for. A fundamental design decision of Fortnet is
not to provide native support for the output files of popular simulation
packages directly. Instead, a separate input format is used and a corresponding
Python class is provided which, based on the Atomic Simulation Environment
(`ASE <https://wiki.fysik.dtu.dk/ase/>`_) that is also implemented in Python,
enables a dataset to be generated easily. To see how you get from the output
files of your simulation package of choice to a Fortnet compatible dataset,
please consult the :ref:`Generating a Dataset <sec-fnetdata>` section.

In this case a file named `training_data` is present in the same folder as the
``fortnet_in.hsd`` input::

  20
  ./dataset/point_01
  ./dataset/point_02
  ./dataset/point_03
       .
       .
       .

The first line contains an integer that specifies the number of ``fnetdata.xml``
paths the current file contains. Following that, the relative (or absolute)
paths to the directories containing the ``fnetdata.xml`` files get listed. Note
that there is no '/' at the end of each path because Fortnet will append the
`/fnetdata.xml` for you. Analogous to the ``Mapping`` block there is an option
(``Standardization``) to perform a simple z-score standardization on the target
values.

In addition, the ``Data`` block also handles the naming scheme of the files
containing all the properties of a single sub-nn of the BPNN, called `netstat`
files in the Fortnet world. The most convenient method, especially for datasets
with multiple atom types, is to use the `Type2FileNames` option. In this case
the only necessary entries are the pre- and suffix of the files and wether to
use lower case characters only (optional, default: No). The parser will then
build appropriate filenames (`./Si.net`, `./C.net`, ...) based on the atom types
found in the dataset at hand. Although not recommended, the output paths and
filenames can also be specified manually, i.e. if different folders are
desired::

  NetstatFiles {
    Si = '/home/user/Silicon.net'
  }


Options
-------
::

  Options {
    Mode = 'train'
    ReadNetStats = No
    RandomSeed = 123456
  }

The basic program behavior gets defined in the ``Option`` block of the input,
starting with the running mode of Fortnet. There are three valid options:
`train`, `validate`, `predict`. As in this example, the `train` mode will
optimize the network with respect to the targets provided by the dataset. A
resumption of the training process based on existing `netstat` files would be
requested by setting the ``ReadNetStats`` entry to `Yes`. To validate the
resulting networks or to predict structures with unknown properties, the
other two modes are used and explained in the
:ref:`First Predictions with Fortnet <sec-firstpredict>` section.

The reproducibility of results is particularly important in scientific fields of
application. To meet this requirement, Fortnet provides a ``RandomSeed`` entry.
By setting a seed you define the initial state of the luxury random number
generator :cite:`ranlux1,ranlux2,ranlux3` that is working in the background and
is responsible for the outcome of the initialization of the sub-nn's and
therefore the training process in general.

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
  |  Fortnet - A BPNN Implementation, Version 0.1                                |
  |                                                                              |
  |  Copyright (C) 2020 - 2021  T. W. van der Heide                              |
  |==============================================================================|

  date: 27.03.2021
  time: 20:57:03, +0100

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

You will also see a list of informationen from the HSD input, as printed below::

  Initialisation

  running in training mode
  random seed: 123456
  read initial netstats: F

  --------------------------------------------------------------------------------

  Sub-NN Details

  hidden layers: 2 2
  activation: tanh

  --------------------------------------------------------------------------------

  ACSF Mappings

  cutoff: 4.0000 Angstrom
  nr. of radial functions: 5
  nr. of angular functions: 4
  Standardization: T

  --------------------------------------------------------------------------------

  Dataset Information

  found: 25 geometries
  in pathfile: training_data
  targets per parameter: .6410

  --------------------------------------------------------------------------------

The entry ``targets per parameter`` is of particular importance. Based on this
ratio you can roughly deduce whether the selected network size is suitable
regarding the dataset that was provided.

Up to this stage of binary execution, the input was parsed and the dataset read.
The ``Calculating ACSF`` statement tells us, that Fortnet has started to map the
structure information to input-suitable ACSF values. As soon as the word `done`
appears, this process is complete and the training process starts::

  Calculating ACSF...done
  Starting training...

       iTrain               Loss          Gradients
  --------------------------------------------------------------------
	1000        0.318759E-02       0.466689E-03
	2000        0.237408E-02       0.397025E-03
	3000        0.196746E-02       0.908079E-04
	4000        0.148653E-02       0.260620E-03
	5000        0.115861E-02       0.250412E-03
  --------------------------------------------------------------------

  Training finished (max. Iterations reached)

While the training process is running, the course of the loss function and the
total gradient of the network parameters are printed regularly, depending on the
``NPrintout`` setting of the ``Training`` block.

Output Files
------------
Depending on the setting of the program behavior in the input file (i.e. running
mode), different output files are created. Running the current example there
will be two files written to disk, appart from the redirected standard output:
`acsf.out`, `Si.net`. The average user does not have to look into either of
these files. They only contain information about the ACSF mappings and the
status of the silicon network, which are necessary for a later resumption of the
training process or for predictions based on the resulting network potential.

In fact, the relevant output ``fnetout.xml`` is only created in validation or
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
