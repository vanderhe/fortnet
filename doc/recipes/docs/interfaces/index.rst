.. _sec-interfaces:

###########################
Interfaces with other codes
###########################



In this section some of the methods for communication between Fortnet and
external software are discussed. This can be used to improve ease of use and
allows you to invoke Fortnet via software that you may be more more familiar
with, such as the Atomic Simulation Environment
(`ASE <https://wiki.fysik.dtu.dk/ase/>`_).

Use of interface communication also offers the possibility of expanding the
applications and functionality of Fortnet. For example, Fortnet can serve as an
energy/force engine for an external driver which then could perform calculations
like geometry optimisation or advanced molecular dynamics based on previously
constructed network potentials.

.. toctree::
   :maxdepth: 1

   sockets/index.rst
   ase/index.rst
