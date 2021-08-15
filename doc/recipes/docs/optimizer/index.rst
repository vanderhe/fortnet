.. _sec-optimizer:
.. highlight:: none

#########
Optimizer
#########

To successively optimize the weight and bias network parameters during the
training iterations, Fortnet provides different algorithms. Depending on the
problem and dataset, the choice of optimizer can have a major impact on
convergence and overall behavior during training. Currently, the following
choices are available:

* Steepest Descent (SD) :cite:`sd`
* Conjugate Gradient (CG) :cite:`cg`
* FIRE (FIRE) :cite:`fire`
* Limited-Memory Broyden–Fletcher–Goldfarb–Shanno (LBFGS) :cite:`lbfgs`


General Optimizer Settings
==========================

Some parameters of the ``Training`` block are universally valid across all
optimizers listed above. The table below lists these entries:

.. list-table::
   :widths: 25 25 25 50
   :header-rows: 0

   * - **Setting**
     - **Type**
     - **Default**
     - **Note**
   * - NIterations
     - Integer
     - Huge()
     - Max. number of training iterations
   * - Threshold
     - Float
     - Tiny()
     - Gradient termination criterion
   * - NPrintout
     - Integer
     - 10
     - Standard output print interval
   * - NSaveNet
     - Integer
     - 100
     - Netstat output save interval
   * - MinDisplacement
     - Float
     - 1e-06
     - Min. displacement in parameters
   * - MaxDisplacement
     - Float
     - 1e+04
     - Max. displacement in parameters
   * - Shuffle
     - Logical
     - No
     - Randomly shuffle order of gradient calculations

Optimizer Specific Settings
===========================

In addition to the universal parameters, there are also optimizer-specific
options. These are the subject of the following sections.

Steepest Descent
----------------

Exemplary HSD ``Training`` block of the ``fortnet_in.hsd`` user input::

  Training = SD {
    Threshold = 1e-08
    NIterations = 10000
    NPrintout = 10
    NSaveNet = 100
    MinDisplacement = 1e-10
    MaxDisplacement = 5e-02
    LearningRate = 0.01
    Shuffle = No
  }

Optimizer specific settings:

.. list-table::
   :widths: 25 25 25 50
   :header-rows: 0

   * - **Setting**
     - **Type**
     - **Default**
     - **Note**
   * - LearningRate
     - Float
     - 0.01
     - uniform weight of gradient components

Conjugate Gradient
------------------

Exemplary HSD ``Training`` block of the ``fortnet_in.hsd`` user input::

  Training = CG {
    Threshold = 1e-08
    NIterations = 10000
    NPrintout = 10
    NSaveNet = 100
    MinDisplacement = 1e-10
    MaxDisplacement = 5e-02
    Shuffle = No
  }

Currently, there are no specific parameters for the conjugate gradient method.

FIRE
----

Exemplary HSD ``Training`` block of the ``fortnet_in.hsd`` user input::

  Training = FIRE {
    Threshold = 1e-08
    NIterations = 10000
    NPrintout = 10
    NSaveNet = 100
    MinDisplacement = 1e-10
    MaxDisplacement = 5e-02
    Shuffle = No
  }

Currently, there are no specific parameters for the conjugate gradient method.

L-BFGS
------

Exemplary HSD ``Training`` block of the ``fortnet_in.hsd`` user input::

  Training = LBFGS {
    Threshold = 1e-08
    NIterations = 10000
    NPrintout = 10
    NSaveNet = 100
    MinDisplacement = 1e-10
    MaxDisplacement = 5e-02
    MaxForQNDisplacement = No
    LineMin = Yes
    Memory = 1000
    Shuffle = No
  }

Optimizer specific settings:

.. list-table::
   :widths: 25 25 25 50
   :header-rows: 0

   * - **Setting**
     - **Type**
     - **Default**
     - **Note**
   * - MaxForQNDisplacement
     - Logical
     - False
     - Consider max. step for quasi-Newton direction
   * - Linemin
     - Logical
     - True
     - Use a line search
   * - Memory
     - Integer
     - 1000
     - Nr. of past iterations to save
