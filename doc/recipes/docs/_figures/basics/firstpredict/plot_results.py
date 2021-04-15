#!/usr/bin/env python3


'''
Script to plot the results of the "First Predictions with Fortnet" section.
'''


import numpy as np
from matplotlib import pyplot as plt


def main():
    '''Main driver routine.'''

    data = np.loadtxt('results.dat')
    nndists = np.arange(2.10, 3.30 + 0.05, 0.05)

    # convert total energy to energy per atom
    data = data / 2.0

    plt.figure(1, figsize=(7, 5))
    plt.title('Comparison of Neural Network Prediction with Targets')
    plt.xlabel(r'Nearest Neighbour Distance [$\mathrm{\AA}$]')
    plt.ylabel('Total Energy [eV] / Atom')

    plt.plot(nndists, data[:, 0], color='blue', label=r'NN',
             linestyle='-', marker='o', markersize=4, zorder=1)
    plt.scatter(nndists, data[:, 1], color='black',
                s=7.0, label=r'DFT', zorder=2)

    plt.tight_layout()
    plt.legend()
    plt.savefig('energy-volume-scan.svg', dpi=900, format='svg')


    data = np.loadtxt('results_plus_extrap.dat')
    nndists = np.arange(2.05, 3.50 + 0.05, 0.05)

    # convert total energy to energy per atom
    data = data / 2.0

    plt.figure(2, figsize=(7, 5))
    plt.title('Comparison of Neural Network Prediction with Targets')
    plt.xlabel(r'Nearest Neighbour Distance [$\mathrm{\AA}$]')
    plt.ylabel('Total Energy [eV] / Atom')

    plt.plot(nndists, data[1:, 0], color='blue', label=r'NN',
             linestyle='-', marker='o', markersize=4, zorder=1)
    plt.scatter(nndists, data[1:, 1], color='black',
                s=7.0, label=r'DFT', zorder=2)
    plt.vlines(2.1, min(data[:, 0]), data[2, 0],
               colors='black', linestyles='dashed')
    plt.vlines(3.3, min(data[:, 0]), data[-5, 0],
               colors='black', linestyles='dashed')

    plt.tight_layout()
    plt.legend()
    plt.savefig('e-v-scan_plus_extrap.svg', dpi=900, format='svg')
    plt.show()


if __name__ == '__main__':
    main()
