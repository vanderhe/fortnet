#!/usr/bin/env python3


'''
Application example of the Fnetdata class, based on a dataset
that provides global system properties as target values to fit.
'''


import os
import numpy as np
from fortformat import Fnetdata
from ase.io.vasp import read_vasp, read_vasp_out


def main():
    '''Main driver routine.'''

    nndists = np.arange(2.00, 3.50 + 0.05, 0.05)

    inpaths = [os.path.join(os.getcwd(), 'vaspdata', entry)
               for entry in sorted(os.listdir('vaspdata'))]

    strucs = []
    energies = np.empty((len(inpaths), 1))

    for ii, inpath in enumerate(inpaths):
        strucs.append(read_vasp(os.path.join(inpath, 'POSCAR')))
        props = read_vasp_out(os.path.join(inpath, 'OUTCAR'))
        energies[ii, 0] = props.get_total_energy()

    fnetdata = Fnetdata(atoms=strucs, globaltargets=energies)
    fnetdata.dump('fnetdata.hdf5')


if __name__ == '__main__':
    main()
