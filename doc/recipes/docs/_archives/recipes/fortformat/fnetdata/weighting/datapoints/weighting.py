#!/usr/bin/env python3


'''
Application example of the Fnetdata class, based on a dataset
that provides global system properties as target values to fit.
The entire dataset is written to a contiguous fnetdata.hdf5 file,
whereby the individual datapoints are weighted by user-input.
'''


import os
import numpy as np
from fortformat import Fnetdata
from ase.io.vasp import read_vasp, read_vasp_out


def main():
    '''Main driver routine.'''

    inpaths = [os.path.join(os.getcwd(), '../../globalTargets/vaspdata', entry)
               for entry in sorted(os.listdir('../../globalTargets/vaspdata'))]

    # start with homogeneous weighting
    weights = np.ones((31,), dtype=int)
    # possibly, certain datapoints are more important
    weights[4:13] = 3

    strucs = []
    energies = np.empty((len(inpaths), 1))

    for ii, inpath in enumerate(inpaths):
        strucs.append(read_vasp(os.path.join(inpath, 'POSCAR')))
        props = read_vasp_out(os.path.join(inpath, 'OUTCAR'))
        energies[ii, 0] = props.get_total_energy()

    fnetdata = Fnetdata(atoms=strucs, targets=energies, atomic=False)
    fnetdata.weights = weights
    fnetdata.dump('fnetdata.hdf5')


if __name__ == '__main__':
    main()
