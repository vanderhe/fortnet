#!/usr/bin/env python3


'''
Application example of a single-point calculation by invoking
Fortnet via the file-IO interface to the ASE Python package.
'''


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
