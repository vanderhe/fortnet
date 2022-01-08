.. highlight:: none
.. _sec-interfaces-ase-fileio:

*******
File-IO
*******

For calculations without heavy file-IO, i.e. systems whose wallclock time is
dominated by the actual network evaluation (e.g. due to large system sizes), the
communication between Fortnet and external software via file-IO may be suitable.

.. note::

   For simulations of small geometries it might be preferable to use a serial
   version of Fortnet, since the overhead in MPI process generation and data
   broadcasting will reverse the speed advantage of an MPI-parallelized build.

Calling Fortnet via ASE
=======================

[Input: `recipes/interfaces/ase/fileio/single/`]

In order for ASE to find the Fortnet executable, an environment variable must be
set. Assuming the use of the BASH shell, this is done as follows (in general an
adaptation to your specific computing environment will be necessary)::

  export FORTNET_COMMAND='~/fortnet/_build/_install/bin/fnet'

An MPI-parallelized build might likewise be configured via::

  export FORTNET_COMMAND='mpirun -np 2 ~/fortnet/_build/_install/bin/fnet'

In this case the single point calculation of a simple :math:`\mathrm{H}_2`
molecule serves as an example. Following the necessary imports, the main script
loads an optimized :math:`\mathrm{H}_2` molecule as provided by the G2-database.
The main method then first instantiates a ``Fortnet()`` calculator object with
options that are crucial for the actual calculation. The prediction run is based
on the previously constructed and converged network potential stored in the
`fortnet.hdf5` file. Subsequently, the calculator gets attached to the system or
geometry respectively and the calculation is started.

Finally, the convergent geometry as well as the forces and corresponding energy
of the system can be read out directly by ASE:

.. code-block:: python

  from ase.build import molecule
  from fnetase import Fortnet

  def main():
      '''Main driver routine.'''

      system = molecule('H2')

      calc = Fortnet(label='H2', atoms=system, restart='fortnet.hdf5')

      calc.calculate(atoms=system, properties=('energy', 'forces'))

      energy = system.get_potential_energy()
      forces = system.get_forces()

  if __name__ == '__main__':
      main()

The script above causes ASE to create an input file (`fortnet_in.hsd`) with the
specified options and invokes Fortnet in the corresponding directory.

By default, the calculator utilizes Fortnet's analytical expressions to infer
atomic forces. If for any reason the numerical calculation via central finite
differences is desired, the keyword argument ``finitediffdelta`` can be set,
which specifies the coordinate shift in angstrom. The calculator will switch to
numerical expressions accordingly:

.. code-block:: python

      calc = Fortnet(label='H2', atoms=system, restart='fortnet.hdf5',
                     finitediffdelta=1e-03)

.. warning::

   For simulations based on the computation of atomic forces, it is strongly
   recommended to use :ref:`regularized networks <sec-loss-regularization>`,
   otherwise the system can/will blow up, especially for geometry optimizations
   and molecular dynamics outlined below.

Geometry Optimization by ASE
============================

[Input: `recipes/interfaces/ase/fileio/geoopt/`]

Apart from the invocation of Fortnet via file-IO, the use of Fortnet as an
energy/force engine in conjunction with an external geometry driver is possible.
To do so, again, set the required environment variable (see explanation above)
and consider a perturbed :math:`\mathrm{H}_2` molecule.

Similar to the previous section, the calculator gets instantiated. Subsequently,
the driver of the geometry optimization is specified and where to write the
trajectory and logfile. In this case, ``BFGS()`` is used, representative of all
the drivers provided by ASE:

.. code-block:: python

  from ase.build import molecule
  from ase.optimize import BFGS
  from fnetase import Fortnet

  def main():
      '''Main driver routine.'''

      system = molecule('H2')

      # perturb H2 molecule
      system.positions[0] = system.positions[0] * 1.3
      system.positions[1] = system.positions[1] * 1.3

      calc = Fortnet(label='H2', atoms=system, restart='fortnet.hdf5')

      system.calc = calc

      opt = BFGS(system, trajectory='opt.traj', logfile='opt.log')
      opt.run(fmax=1.0e-006)

      energy = system.get_potential_energy()
      forces = system.get_forces()

  if __name__ == '__main__':
      main()

The script shown causes ASE to generate appropriate input files for each step of
the geometry optimization. Note that this can lead to heavy file-IO and thus a
significant increase in wallclock time, depending on the speed of the storage
used. Therefore it is advisable to perform such calculations on a ramdisk, if
available.

Molecular Dynamics by ASE
=========================

[Input: `recipes/interfaces/ase/fileio/md/`]

Apart from the invocation of Fortnet via file-IO, the use of Fortnet as an
energy/force engine in conjunction with an external molecular dynamics driver is
possible. To do so, again, set the required environment variable (see
explanation above) and consider an :math:`\mathrm{H}_2` molecule from the
G2-database.

Similar to the previous section, the calculator gets instantiated. Subsequently,
a Maxwell-Boltzmann distribution is used to initialize atomic velocities at the
desired temperature. A canonical ensemble (constant NVT) is employed to keep the
volume and temperature fixed during the simulation, representative for all
ensembles provided by ASE. A simple function prints atomic energies at a desired
interval.

.. code-block:: python

  from ase.build import molecule
  from ase.io.trajectory import Trajectory
  from ase.md.langevin import Langevin
  from ase.md.velocitydistribution import (MaxwellBoltzmannDistribution,
					   Stationary, ZeroRotation)
  from ase import units
  from fnetase import Fortnet

  def main():
      '''Main driver routine.'''

      system = molecule('H2')

      system.calc = Fortnet(label='H2', atoms=system, restart='fortnet.hdf5')

      MaxwellBoltzmannDistribution(system, temperature_K=200)
      Stationary(system)
      ZeroRotation(system)

      dyn = Langevin(system, 1.0 * units.fs, friction=1e-02, temperature_K=200)

      def printenergy(atoms=system):
	  '''Prints the potential, kinetic and total energy.'''
	  epot = atoms.get_potential_energy() / len(atoms)
	  ekin = atoms.get_kinetic_energy() / len(atoms)
	  print('Energy per atom: Epot = %.3feV  Ekin = %.3feV (T=%3.0fK)  '
		'Etot = %.3feV' % (epot, ekin, ekin / (1.5 * units.kB),
				   epot + ekin))

      dyn.attach(printenergy, interval=10)

      traj = Trajectory('md.traj', 'w', system)
      dyn.attach(traj.write, interval=1)

      printenergy()
      dyn.run(200)

  if __name__ == '__main__':
      main()

After running the simulation, you can study the result with the command::

  ase gui md.traj
