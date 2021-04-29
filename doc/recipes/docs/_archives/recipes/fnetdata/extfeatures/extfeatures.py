#!/usr/bin/env python3


'''
Application example of the Fortformat class, based on a dataset
that provides global system properties as target values to fit.
The dataset is extended by user specified external atomic features.
'''


import os
import numpy as np
from fortformat import Fortformat
from ase.io.vasp import read_vasp, read_vasp_out


def main():
    '''Main driver routine.'''

    np.random.seed(42)

    inpaths = [os.path.join(os.getcwd(), '../globalTargets/vaspdata', entry)
               for entry in sorted(os.listdir('../globalTargets/vaspdata'))]

    strucs = []
    features = []
    energies = np.empty((len(inpaths), 1))

    for ii, inpath in enumerate(inpaths):
        struc = read_vasp(os.path.join(inpath, 'POSCAR'))
        strucs.append(struc)
        props = read_vasp_out(os.path.join(inpath, 'OUTCAR'))
        energies[ii, 0] = props.get_total_energy()
        features.append(np.random.random_sample((len(struc), 3)))

    fnetdata = Fortformat(strucs, 'fnetdata.xml', targets=energies,
                          features=features, atomic=False, frac=True)
    fnetdata.dump()


if __name__ == '__main__':
    main()
