.. _sec-acsf:
.. highlight:: none

################################
Atom-Centered Symmetry Functions
################################

A fundamental problem with the use of neural networks for the construction of
high-dimensional potential-energy hypersurfaces is the generation of suitable
input features with which the networks can be fed. The raw coordinates are
unsuitable for several reasons, as pointed out by J. Behler :cite:`acsf`:
Translations or rotations of a system are energy-conserving operations,
regardless of the change in the absolute coordinates. However, this condition
wouldn't be strictly fulfilled for networks that are trained on the absolute
atomic coordinates. Otherwise the network would be fitted to potentially
contradicting data, since two different inputs would have to lead to exactly the
same output (i.e. total energy). The same reasoning applies to the exchange of
two atoms of the same type as well as the requirement for a constant number of
input features, regardless of possibly changing coordination numbers
(i.e. during MD simulations).

Therefore Fortnet implements the so-called Atom-centered Symmetry Functions
(ACSF) :cite:`acsf`, proposed by J. Behler as a possible choice of symmetry
mappings, listed below:

.. math::

   \begin{align*}
   f_\mathrm{c}(R_{ij}) &=
   \begin{cases}
   \frac{1}{2}\left[\cos\left(\frac{\pi R_{ij}}{R_\mathrm{c}}\right)+
   1\right]&\text{ for }R_{ij}\leq R_\mathrm{c} \\ 0 &\text{ for }R_{ij}>
   R_\mathrm{c}
   \end{cases}
   \end{align*}

radial functions:

.. math::

   \begin{align*}
   G_i^1 &= \sum_j f_\mathrm{c}(R_{ij}) \\
   G_i^2 &= \sum_j e^{-\eta (R_{ij} - R_\mathrm{s})^2}\cdot f_\mathrm{c}(R_{ij})
   \\
   G_i^3 &= \sum_j \cos(\kappa R_{ij})\cdot f_\mathrm{c}(R_{ij})
   \end{align*}

.. math::

   \Theta_{ijk} =
   \arccos\left(\frac{\boldsymbol{R}_{ij}\cdot
   \boldsymbol{R}_{ik}}{R_{ij}\cdot R_{ik}}\right)

angular functions:

.. math::

   \begin{align*}
   G_i^4 &= 2^{1-\xi} \sum_{j, k\neq i}^{\mathrm{all}}(1+\lambda\cos
   \Theta_{ijk})^\xi \cdot e^{-\eta (R_{ij}^2 + R_{ik}^2 + R_{jk}^2)}\cdot
   f_\mathrm{c}(R_{ij})f_\mathrm{c}(R_{ik})f_\mathrm{c}(R_{jk}) \\
   G_i^5 &= 2^{1-\xi} \sum_{j, k\neq i}^{\mathrm{all}} (1 + \lambda\cos
   \Theta_{ijk})^\xi \cdot e^{-\eta (R_{ij}^2 + R_{ik}^2)}\cdot
   f_\mathrm{c}(R_{ij})f_\mathrm{c}(R_{ik})
   \end{align*}


Please note that currently only an automatic generation scheme is available via
the HSD input. An explanation of this parameter generation can be found in the
following section.

==============================
Automatic Parameter Generation
==============================

The automatic ACSF parameter generation scheme aims to cover the cutoff sphere
as evenly as possible, with the number of symmetry functions available
:cite:`prophet`. The user HSD input is already known from the introductory
sections and looks similar to this::

  Mapping = ACSF {
    NRadial = 10
    NAngular = 8
    RCut = 4.0
  }

The corresponding parameters are then generated based on the number ``NRadial``
of :math:`G^2` and the number ``NAngular`` of :math:`G^5` functions specified,
as well as the cutoff radius ``RCut`` in Angstrom.

Radial Parameters
-----------------

The :math:`G^2` function is used as radial mapping, therefore sensible values
for the peak position :math:`R_\mathrm{s}` and width :math:`\eta` of the
Gaussians are to be determined. :math:`R_\mathrm{s}` is chosen to go from 0 to
:math:`R_\mathrm{c}` in equidistant steps. The width evaluates to the following
more or less arbitrary equation:

.. math::

   \eta = \frac{5\cdot \ln(10)}{4a^2}\quad\text{where}\quad
   a = \frac{R_\mathrm{c}}{N_\mathrm{rad}-1}

The figure below shows the :math:`G^2` functions resulting for the HSD input
example above:

.. figure:: ../_figures/acsf/g2.svg
   :width: 100%
   :align: center
   :alt: Exemplary Visualization of :math:`G^2` Functions.



Angular Parameters
------------------

The :math:`G^5` function is used as angular mapping, therefore sensible values
for the pre-factor exponent :math:`\xi`, :math:`\lambda`-parameter and width
:math:`\eta` of the Gaussians are to be determined. :math:`\xi` is chosen to go
from 1 to 16 in equidistant steps, in order to obtain functions that are
strictly separated from one another. The width, again, evaluates to the
following more or less arbitrary equation:

.. math::

   \eta = \frac{2\cdot \ln(10)}{R_\mathrm{c}^2}

Last but not least, the parameter :math:`\lambda` alternately takes the values
-1 and 1, as suggested by the original ACSF paper :cite:`acsf`:

.. math::

   \lambda = \{-1,1\}

=============================
Atom-Specific Scaling Factors
=============================

One possible extension of the ACSF mappings is the incorporation of external
atomic scaling factors :math:`q_j`, hereinafter referred to as atom identifiers,
in the cutoff function :math:`f_\mathrm{c}`:

.. math::

   \begin{align*}
   f_\mathrm{c}(R_{ij}) &=
   \begin{cases}
   \frac{\color{red}{q_j}}{2}\left[\cos\left(\frac{\pi R_{ij}}{R_\mathrm{c}}
   \right)+1\right]&\text{ for }R_{ij}\leq R_\mathrm{c} \\ 0 &\text{ for }R_{ij}
   >R_\mathrm{c}
   \end{cases}
   \end{align*}

In addition to the structural information, atom-specific
features can thus be taken into account. A reasonable choice for the atom
identifiers would be, for example, the Mulliken populations from a quantum
mechanical simulation of the system.

The atom identifiers are extracted from the external atomic features provided by
the dataset. To do so, the dataset must provide at least one external feature
per atom to choose and an appropriate entry in the ``External`` HSD-block of the
user input has to be set::

   External {
         .
	 .
	 .
     AtomID = 1
   }

The integer index corresponds to the atomic feature of the dataset to extract
and must take a value between one and the number of features per atom provided.

=====================
Multi-Species Systems
=====================

A modification of the ACSF can be useful for datasets that contain systems with
several atom types. In the form as listed above, no distinction is made between
different species, which is why the resulting values do not depend on the
species of the atoms contained in the cutoff sphere, defined by
:math:`R_\mathrm{c}`.

To overcome this limitation, Fortnet offers the possibility to specify
so-called species identifiers in the `External` block of the ``fortnet_in.hsd``
input file::

   External {
         .
	 .
	 .
     SpeciesID {
       C = 0.364302
       Si = 0.247609
     }
   }

The default value of the species identifier is 1.0. Element-specific properties
such as the Hubbard :math:`U` of the DFTB method are particularly suitable. The
identifiers are then used as a pre-factor in the cutoff function
:math:`f_\mathrm{c}` and are therefore included in all ACSF:

.. math::

   \begin{align*}
   f_\mathrm{c}(R_{ij}) &=
   \begin{cases}
   \frac{\color{red}{U_j}}{2}\left[\cos\left(\frac{\pi R_{ij}}{R_\mathrm{c}}
   \right)+1\right]&\text{ for }R_{ij}\leq R_\mathrm{c} \\ 0 &\text{ for }R_{ij}
   >R_\mathrm{c}
   \end{cases}
   \end{align*}

These species identifiers are compatible with the atom identifiers. If both are
provided, the corresponding scalings are multiplied in the cutoff function. 
