.. _sec-loss:
.. highlight:: none

##############
Loss Functions
##############

Choosing a suitable loss function for the task at hand can make a significant
difference, primarily depending on the dataset targets. In order to give you a
certain freedom in this sense (i.e. when it comes to weighting outliers),
Fortnet implements the following functions:

  - mean squared loss (mse)
  - root mean square loss (rms)
  - mean absolute loss (mae)
  - mean absolute percentage loss (mape)

The loss function used during the training is selected in the ``Training`` block
of the HSD input. The functions and an associated, exemplary training block, are
listed below, assuming a dataset of :math:`N` targets :math:`y_i^\mathrm{ref}`
and network predictions :math:`y_i^\mathrm{nn}`.

.. note::
   The default loss function is the mean squared error (``mse``). 

Mean Squared Error
==================
.. math::

  \begin{align*}
  C = \frac{1}{N}\sum_{i=1}^N \left(y_i^\mathrm{ref} - y_i^\mathrm{nn}\right)^2
  \end{align*}

::

  Training = LBFGS {
        .
	.
	.
    Loss = 'mse'
  }

Root Mean Square Error
======================
.. math::

  \begin{align*}
  C = \sqrt{\frac{1}{N}\sum_{i=1}^N \Big(y_i^\mathrm{ref} -
  y_i^\mathrm{nn}\Big)^2}
  \end{align*}

::

  Training = LBFGS {
        .
	.
	.
    Loss = 'rms'
  }

Mean Absolute Error
===================
.. math::

  \begin{align*}
  C = \frac{1}{N}\sum_{i=1}^N |y_i^\mathrm{ref} - y_i^\mathrm{nn}|
  \end{align*}

::

  Training = LBFGS {
        .
	.
	.
    Loss = 'mae'
  }

Mean Absolute Percentage Error
==============================
.. math::

  \begin{align*}
  C = \frac{100}{N}\sum_{i=1}^N \frac{|y_i^\mathrm{ref} - y_i^\mathrm{nn}|}
  {|y_i^\mathrm{ref}|}
  \end{align*}

::

  Training = LBFGS {
        .
	.
	.
    Loss = 'mape'
  }
