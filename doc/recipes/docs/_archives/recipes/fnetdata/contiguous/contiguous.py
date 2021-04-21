#!/usr/bin/env python3


'''
Application example of the Fortformat class, based on a dataset
that provides global system properties as target values to fit.
The entire dataset is written to a contiguous fnetdata.xml file.
'''


import os
import numpy as np
from fortformat import Fortformat
from ase.io.vasp import read_vasp, read_vasp_out


def main():
    '''Main driver routine.'''

    nndists = np.arange(2.00, 3.50 + 0.05, 0.05)

    inpaths = [os.path.join(os.getcwd(), 'vaspdata', entry)
               for entry in sorted(os.listdir('vaspdata'))]
    outpaths = [os.path.join(os.getcwd(), 'dataset', 'nndist_{:.3f}'
                             .format(nndist)) for nndist in nndists]

    strucs = []
    energies = np.empty((len(inpaths), 1))

    for ii, inpath in enumerate(inpaths):
        strucs.append(read_vasp(os.path.join(inpath, 'POSCAR')))
        props = read_vasp_out(os.path.join(inpath, 'OUTCAR'))
        energies[ii, 0] = props.get_total_energy()

    fnetdata = Fortformat(strucs, 'fnetdata.xml', targets=energies,
                          atomic=False, frac=True)
    fnetdata.dump()


if __name__ == '__main__':
    main()
