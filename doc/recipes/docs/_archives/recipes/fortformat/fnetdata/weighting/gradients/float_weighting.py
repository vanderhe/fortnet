#!/usr/bin/env python3


'''
Application example of the Fnetdata class, based on a dataset
that provides global system properties as target values to fit.
The entire dataset is written to a contiguous fnetdata.hdf5 file,
whereby atomic gradient contributions are weighted by user-input.
'''


import os
import numpy as np
from fortformat import Fnetdata
from ase.io.vasp import read_vasp, read_vasp_out


def main():
    '''Main driver routine.'''

    inpaths = [os.path.join(os.getcwd(), '../../globalTargets/vaspdata', entry)
               for entry in sorted(os.listdir('../../globalTargets/vaspdata'))]

    # fix random seed for reproduction purposes
    np.random.seed(42)

    strucs = []
    atomicweights = []
    energies = np.empty((len(inpaths), 1))

    for ii, inpath in enumerate(inpaths):
        struc = read_vasp(os.path.join(inpath, 'POSCAR'))
        strucs.append(struc)
        props = read_vasp_out(os.path.join(inpath, 'OUTCAR'))
        energies[ii, 0] = props.get_total_energy()
        # float-valued atomic gradient weighting in interval [1, 10]
        atomicweights.append(np.asfarray(
            np.random.randint(1, 10, len(struc), dtype=int)))

    fnetdata = Fnetdata(atoms=strucs, targets=energies, atomic=False)
    fnetdata.atomicweights = atomicweights
    fnetdata.dump('fnetdata_float_weighting.hdf5')


if __name__ == '__main__':
    main()
