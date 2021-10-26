#!/usr/bin/env python3


'''
Application example of an H2 molecular dynamics by invoking
Fortnet via the file-IO interface to the ASE Python package.
'''


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

    system.calc = Fortnet(label='H2', atoms=system, restart='fortnet.hdf5',
                          finiteDiffDelta=1e-02)

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
