.. _sec-extfeatures:
.. highlight:: none

**************************************
Incorporating External Atomic Features
**************************************

[Input: `recipes/basics/extfeatures/`]

The possible incorporation of user-specified, external atomic features is
arguably one of the most important capabilities of `Fortnet`, since this kind of
software can never claim that all conceivable feature generators are
implemented natively. In order to give the user the greatest possible freedom in
generating additional input features (as a supplement to the
:ref:`ACSF <sec-acsf>`), `Fortformat` and `Fortnet` support external atomic
features. After reading this short section, you will be able to make use of it.

Providing the Input
===================

As an input, we will use conventional HSD blocks, that we are already familiar
with from the previous sections::

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

  Training = SD {
    Threshold = 1e-08
    NIterations = 5000
    NPrintout = 1000
    NSaveNet = 1000
    MinDisplacement = 1e-10
    MaxDisplacement = 5e-02
    Loss = 'rms'
  }

  Data {
    Dataset = 'fnetdata.xml'
    Standardization = No
    NetstatFiles = Type2FileNames {
      Prefix = "./"
      Suffix = ".net"
      LowerCaseTypeName = No
    }
  }

  External {
    Features = 1 2 3
  }

  Options {
    Mode = 'train'
    ReadNetStats = No
    RandomSeed = 123456
  }

The only addition is located in the ``External`` block. By specifying the
``Features`` entry we ask Fortnet to incorporate the corresponding entries of
the dataset, provided they are available. So let's take a closer look at the
specification options in the input. In this case, a dataset with three external
features per atom is used. A dedicated :ref:`section <sec-fnetdata_extFeatures>`
introduces the generation of such a dataset, using the ``Fortformat`` Python
class.

There are essentially two input formats for specifying the external features
to be used:

1. Explicit specification of the feature indices as they appear in the dataset::

       External {
          .
	  .
	 Features = 1 2 3
       }

  The advantage is that certain features can be left out or even used multiple
  times like this::

       External {
          .
	  .
	 Features = 1 1 2 2 2
       }

2. Specification of a contiguous range of indices (bounds included)::

       External {
          .
	  .
	 Features = 1:3
       }

  This format can therefore be used to specify all the entries in the dataset in
  a convenient way::

       External {
          .
	  .
	 Features = 1:-1
       }

  Or for example up to the penultimate entry: ``Features = 1:-2``


Standard Output
===============

As far as the standard output is concerned, an additional entry will appear on
the screen, which provides information about the number and indexing of the
features used::

  --------------------------------------------------------------------------------

  External Features

  nr. of external features: 3
  dataset indices: 1 2 3

  --------------------------------------------------------------------------------
