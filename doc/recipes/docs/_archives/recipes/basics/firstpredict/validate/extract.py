#!/usr/bin/env python3


'''
Application example of the Fortformat package, based on an output
file that contains network predictions and corresponding targets.
'''


import numpy as np
import matplotlib.pyplot as plt
from fortformat import Fnetout


def main():
    '''Main driver routine.'''

    nndists = np.arange(2.10, 3.30 + 0.05, 0.05)

    fnetout = Fnetout('fnetout.hdf5')
    predictions = fnetout.predictions
    targets = fnetout.targets

    plt.figure(figsize=(7, 5))
    plt.title('Comparison of Neural Network Predictions with Targets')
    plt.xlabel(r'Nearest Neighbour Distance [$\mathrm{\AA}$]')
    plt.ylabel('Total Energy [eV / Atom]')

    plt.plot(nndists, predictions / 2.0, color='blue', label='NN')
    plt.scatter(nndists, targets / 2.0, s=10, color='black', label='DFT')

    plt.tight_layout()
    plt.legend()
    plt.savefig('comparison.svg', dpi=900, format='svg')


if __name__ == '__main__':
    main()
