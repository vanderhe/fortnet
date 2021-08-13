#!/usr/bin/env python3
'''
Create fnetdata.xml files out of existing dataset.

@author: Tammo van der Heide
'''


import os
import numpy as np
from fortformat import Fortformat
from ase.io.vasp import read_vasp


def get_free_energy(filename):
    ''''''

    keyword = 'free  energy   TOTEN  ='
    with open(filename, 'r') as outfile:
        for line in outfile:
            if keyword in line:
                content = line

    return float(content.strip().strip(keyword).split()[0])


def main():
    '''Main driver routine.'''

    with open('inpaths', 'r') as infile:
        inpaths = infile.readlines()

    inpaths = [entry.strip('\n') for entry in inpaths]

    atoms = []
    energies = np.empty((len(inpaths), 1))
    weights = np.empty((len(inpaths)), dtype=int)

    for ii, inpath in enumerate(inpaths):
        struc = read_vasp(os.path.join(inpath, 'POSCAR'))
        atoms.append(struc)

        weights[ii] = ii + 1
        energies[ii, 0] = get_free_energy(os.path.join(inpath, 'OUTCAR'))

    fnetdata = Fortformat(atoms=atoms, targets=energies, atomic=False)
    fnetdata.weights = weights
    fnetdata.dump('fnetdata.hdf5')


main()
