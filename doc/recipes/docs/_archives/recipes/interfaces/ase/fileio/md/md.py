#!/usr/bin/env python3


'''
Application example of an Si64 molecular dynamics by invoking
Fortnet via the file-IO interface to the ASE Python package.
'''


from ase.io import read
from ase.io.trajectory import Trajectory
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase import units
from fnetase import Fortnet


def main():
    '''Main driver routine.'''

    system = read('POSCAR')

    system.calc = Fortnet(label='Si64', atoms=system, restart='fortnet.hdf5',
                          finiteDiffDelta=1e-02)

    MaxwellBoltzmannDistribution(system, temperature_K=300)

    dyn = Langevin(system, 1.0 * units.fs, friction=1e-02, temperature_K=300)

    def printenergy(atoms=system):
        '''Prints the potential, kinetic and total energy.'''
        epot = atoms.get_potential_energy() / len(atoms)
        ekin = atoms.get_kinetic_energy() / len(atoms)
        print('Energy per atom: Epot = %.3feV  Ekin = %.3feV (T=%3.0fK)  '
              'Etot = %.3feV' % (epot, ekin, ekin / (1.5 * units.kB),
                                 epot + ekin))

    dyn.attach(printenergy, interval=1)

    traj = Trajectory('md.traj', 'w', system)
    dyn.attach(traj.write, interval=1)

    printenergy()
    dyn.run(10)


if __name__ == '__main__':
    main()
