#!/usr/bin/env python3


'''
Application example of an H2 geometry optimization by invoking
Fortnet via the file-IO interface to the ASE Python package.
'''


from ase.build import molecule
from ase.optimize import BFGS
from fnetase import Fortnet


def main():
    '''Main driver routine.'''

    system = molecule('H2')

    # perturb H2 molecule
    system.positions[0] = system.positions[0] * 1.3
    system.positions[1] = system.positions[1] * 1.3

    calc = Fortnet(label='H2', atoms=system, restart='fortnet.hdf5',
                   finiteDiffDelta=1e-03)

    system.calc = calc

    opt = BFGS(system, trajectory='opt.traj', logfile='opt.log')
    opt.run(fmax=1.0e-006)

    energy = system.get_potential_energy()
    print(energy)
    forces = system.get_forces()
    print(forces)


if __name__ == '__main__':
    main()
