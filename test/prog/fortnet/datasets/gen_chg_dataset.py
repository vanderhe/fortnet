#!/usr/bin/env python3
'''
Create fnetdata.xml files out of existing dataset.

@author: Tammo van der Heide
'''


import os
import random
import numpy as np
from fortformat import Fortformat
from ase.io.vasp import read_vasp


Bohr__AA = 0.529177249
AA__Bohr = 1.0 / Bohr__AA
Hartree__eV = 27.2113845
eV__Hartree = 1.0 / Hartree__eV


def get_repulsive_energy(filename):
    '''Extracts the repulsive energy out of detailed.out'''

    keyword = 'Repulsive energy:'
    with open(filename, 'r') as detailed_out:
        for line in detailed_out:
            if keyword in line:
                elen_line = line

    elen = elen_line.strip().strip(keyword).split()[0]

    return float(elen) * Hartree__eV


def get_mulliken_population(filename, natoms):
    '''Extracts the atomic Mulliken populations out of detailed.out'''

    population = np.empty((natoms, 3))

    iline = 0

    keyword = 'Atom populations'
    with open(filename, 'r') as detailed:
        for line in detailed:
            iline = iline + 1
            if keyword in line:
                break

    population[:, 0] = np.loadtxt(filename, skiprows=iline+1, max_rows=natoms)[:, 1]

    iline = 0

    keyword = 'l-shell populations'
    with open(filename, 'r') as detailed:
        for line in detailed:
            iline = iline + 1
            if keyword in line:
                break

    population[:, 1] = np.loadtxt(filename, skiprows=iline+1, max_rows=2*natoms)[0::2, 3]
    population[:, 2] = np.loadtxt(filename, skiprows=iline+1, max_rows=2*natoms)[1::2, 3]

    return population


def main():
    '''Main driver routine.'''

    with open('inpaths', 'r') as infile:
        inpaths = infile.readlines()

    inpaths = [entry.strip('\n') for entry in inpaths]

    strucs = []
    atomic_nums = []
    energies = np.empty((len(inpaths), 1))
    weights = np.empty((len(inpaths)), dtype=int)

    for ii, inpath in enumerate(inpaths):
        struc = read_vasp(os.path.join(inpath, 'POSCAR'))
        strucs.append(struc)

        repen = get_repulsive_energy(os.path.join(inpath, 'detailed.out'))
        zz = np.empty((len(struc), 1), dtype=float)
        zz[:, 0] = struc.get_atomic_numbers()
        atomic_nums.append(zz)

        weights[ii] = ii + 1
        energies[ii, 0] = repen

    fnetdata = Fortformat(atoms=strucs, targets=energies,
                          atomic=False)
    fnetdata.weights = weights
    fnetdata.dump('fnetdata.hdf5')


main()


# def main():
#     '''Main driver routine.'''

#     with open('inpaths_part_63', 'r') as infile:
#         inpaths_63 = infile.readlines()

#     with open('inpaths_part_64', 'r') as infile:
#         inpaths_64 = infile.readlines()

#     inpaths_63 = [entry.strip('\n') for entry in inpaths_63]
#     inpaths_64 = [entry.strip('\n') for entry in inpaths_64]

#     strucs = []
#     populations_63 = []
#     populations_64 = []
#     energies = np.empty((len(inpaths_63) + len(inpaths_64), 1))

#     for ii, inpath in enumerate(inpaths_64):
#         struc = read_vasp(os.path.join(inpath, 'POSCAR'))
#         natoms = len(struc)
#         strucs.append(struc)

#         population = get_mulliken_population(
#             os.path.join(inpath, 'detailed.out'), natoms)
#         # population = np.empty((natoms, 3), dtype=float)
#         # population[:, 1] = np.random.random_sample((natoms,)) * 4.0
#         # population[:, 2] = 4.0 - population[:, 1]
#         # population[:, 0] = population[:, 1] + population[:, 2]
#         populations_64.append(population)
#         repen = get_repulsive_energy(os.path.join(inpath, 'detailed.out'))

#         energies[ii, 0] = repen

#     for ii, inpath in enumerate(inpaths_63):
#         struc = read_vasp(os.path.join(inpath, 'POSCAR'))
#         natoms = len(struc)
#         strucs.append(struc)

#         population = get_mulliken_population(
#             os.path.join(inpath, 'detailed.out'), natoms)
#         populations_63.append(population)
#         repen = get_repulsive_energy(os.path.join(inpath, 'detailed.out'))

#         energies[ii + len(inpaths_64), 0] = repen

#     for jj in range(10):
#         random.shuffle(populations_63)
#         random.shuffle(populations_64)

#     populations = []
#     populations += populations_64
#     populations += populations_63

#     fnetdata = Fortformat(strucs, 'fnetdata.xml', targets=energies,
#                           features=populations, atomic=False, frac=True)
#     fnetdata.dump()


# main()
