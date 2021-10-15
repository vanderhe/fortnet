.. _sec-transfer:
.. highlight:: none

####################
Activation Functions
####################

Choosing a suitable activation function for the task at hand can make a
significant difference, depending on the situation or dataset. In order to give
you a certain freedom in this sense, Fortnet implements the following functions:

  - hyperbolic tangent
  - arcus tangent
  - sigmoid function
  - softplus function
  - gaussian function
  - (leaky) ReLU function
  - Bent identity function
  - heaviside function
  - linear function

The activation function is selected in the ``Network`` block of the HSD input.
The functions and an associated, exemplary network block, are listed below.

.. note::
   A linear function is always used as the activation function of the output
   layer, in order to be able to represent all real-valued results.


Hyperbolic Tangent
==================
::

  Network = BPNN {
    Hidden = 2 2
    Activation = tanh
  }

.. figure:: ../_figures/transfer/tanh.svg
   :width: 100%
   :align: center
   :alt: Plot of the hyperbolic tangent.

Arcus Tangent
==================
::

  Network = BPNN {
    Hidden = 2 2
    Activation = atan
  }

.. figure:: ../_figures/transfer/atan.svg
   :width: 100%
   :align: center
   :alt: Plot of the arcus tangent.

Sigmoid
=======
::

  Network = BPNN {
    Hidden = 2 2
    Activation = sigmoid
  }

.. figure:: ../_figures/transfer/sigmoid.svg
   :width: 100%
   :align: center
   :alt: Plot of sigmoid activation function.

SoftPlus
========
::

  Network = BPNN {
    Hidden = 2 2
    Activation = softplus
  }

.. figure:: ../_figures/transfer/softplus.svg
   :width: 100%
   :align: center
   :alt: Plot of softplus activation function.

Gaussian
========
::

  Network = BPNN {
    Hidden = 2 2
    Activation = gaussian
  }

.. figure:: ../_figures/transfer/gaussian.svg
   :width: 100%
   :align: center
   :alt: Plot of gaussian activation function.

ReLU
====
::

  Network = BPNN {
    Hidden = 2 2
    Activation = relu
  }

.. figure:: ../_figures/transfer/relu.svg
   :width: 100%
   :align: center
   :alt: Plot of relu activation function.

Leaky ReLU
==========
::

  Network = BPNN {
    Hidden = 2 2
    Activation = lrelu
  }

.. figure:: ../_figures/transfer/lrelu.svg
   :width: 100%
   :align: center
   :alt: Plot of leaky ReLU activation function.

Bent Identity
=============
::

  Network = BPNN {
    Hidden = 2 2
    Activation = bent
  }

.. figure:: ../_figures/transfer/bent.svg
   :width: 100%
   :align: center
   :alt: Plot of Bent identity activation function.

Heaviside
=========
::

  Network = BPNN {
    Hidden = 2 2
    Activation = heaviside
  }

.. figure:: ../_figures/transfer/heaviside.svg
   :width: 100%
   :align: center
   :alt: Plot of heaviside activation function.

Linear
======
::

  Network = BPNN {
    Hidden = 2 2
    Activation = linear
  }

.. figure:: ../_figures/transfer/linear.svg
   :width: 100%
   :align: center
   :alt: Plot of linear activation function.
