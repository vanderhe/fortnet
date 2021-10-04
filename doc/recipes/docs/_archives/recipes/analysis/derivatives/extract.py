#!/usr/bin/env python3

'''
Application example of the Fortformat package, based on an output
file that contains atomic forces, resulting from finite differences.
'''


import numpy as np
from fortformat import Fnetout


def main():
    '''Main driver routine.'''

    fnetout = Fnetout('fnetout.hdf5')
    forces = fnetout.forces

    # print forces of each datapoint, network
    # output and atom to illustrate the indexing:
    for idata in range(len(forces)):
        for iout in range(len(forces[idata])):
            for iatom in range(np.shape(forces[idata][iout])[0]):
                print(forces[idata][iout][iatom])


if __name__ == '__main__':
    main()
