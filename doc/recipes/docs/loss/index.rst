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
    Loss = mse
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
    Loss = rms
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
    Loss = mae
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
    Loss = mape
  }

Regularization
==============

An equally simple and effective method to prevent overfitting while training a
neural network is loss-based regularization. By adding an additional penalty
term :math:`\tilde{C}` to the base loss :math:`C_0`, the assembly of a spiky
hypersurface, due to high connection weights, can be mitigated:

.. math::

  \begin{align*}
  C = C_0 + \lambda\tilde{C}
  \end{align*}

The ``Strength`` of the penalty is regulated by the :math:`\lambda`-parameter.
Fortnet supports :math:`L_1` (lasso) and :math:`L_2` (ridge) regularizations as
well as a mixture of both (elastic net), serving different purposes. The
specification takes place in the ``Training`` block of the HSD input, e.g.::

  Training = LBFGS {
        .
	.
	.
      Regularization = Ridge {
        Strength = 1.0
      }
  }

More detailed descriptions of each variant are given below.

L1 - Lasso
----------
Lasso regression adds a penalty based on the raw magnitudes of weight
coefficients. It is often referred to as :math:`L_1` regularization and might
lead to a certain feature selection by zeroing out some of the weights:

.. math::

  \begin{align*}
  \tilde{C} &= \frac{1}{n}\sum_{i=1}^{n}|w_i| \\
  \frac{\partial\tilde{C}}{\partial w_i} &= \frac{\mathrm{sgn}(w_i)}{n}
  \end{align*}

The specification takes place in the ``Training`` block of the HSD input, e.g.::

  Training = LBFGS {
        .
	.
	.
      Regularization = Lasso {
        Strength = 1.0
      }
  }

L2 - Ridge
----------
Ridge regression adds a penalty on particularly large weight coefficients. It is
often referred to as :math:`L_2` regularization and, analogous to lasso
regression, shrinks the weights and reduces model complexity:

.. math::

  \begin{align*}
  \tilde{C} &= \frac{1}{2n}\sum_{i=1}^{n}w_i^2 \\
  \frac{\partial\tilde{C}}{\partial w_i} &= \frac{w_i}{n}
  \end{align*}

The specification takes place in the ``Training`` block of the HSD input, e.g.::

  Training = LBFGS {
        .
	.
	.
      Regularization = Ridge {
        Strength = 1.0
      }
  }

Elastic Net
-----------
Elastic net regularization mixes :math:`L_1` and :math:`L_2` contributions in a
certain ratio, determined by the :math:`\alpha`-parameter:

.. math::

  \begin{align*}
  \tilde{C} &= \frac{1}{n}\left(\frac{1-\alpha}{2}\cdot\sum_{i=1}^{n}w_i^2
  + \alpha\cdot\sum_{i=1}^{n}|w_i|\right) \\
  \frac{\partial\tilde{C}}{\partial w_i} &= \frac{1}{n}\left[(1-\alpha)\cdot w_i
  + \alpha\cdot\mathrm{sgn}(w_i)\right]
  \end{align*}

For :math:`\alpha = 0` or :math:`\alpha = 1` this results in :math:`L_2` or
:math:`L_1` regularization respectively. The specification takes place in the
``Training`` block of the HSD input, e.g.::

  Training = LBFGS {
        .
	.
	.
    Regularization = ElasticNet {
      Strength = 1.0
      Alpha = 0.5
    }
  }
