#!/usr/bin/env python3


'''
Visualization of the radial Gaussian G2 functions.

The example contains 10 radial functions,
limited by a cutoff of 4.0 angstroms.
'''


import numpy as np
from matplotlib import pyplot as plt


RCUT = 4.0
NRADIAL = 10


def main():
    '''Main driver routine.'''


    rr = np.linspace(0.0, RCUT, 1000)
    rs = np.linspace(0.0, RCUT, NRADIAL)
    aa = RCUT / (NRADIAL - 1.0)
    eta = 5.0 * np.log(10.0) / (4.0 * aa**2)

    plot_g2(rr, rs, eta)


def plot_g2(rr, rs, eta):
    '''Plots radial Gaussian G2 functions for given parameters.

    Args:

        rr (1darray): distances between atoms
        rs (1darray): Rs parameters of G2 functions
        eta (float): eta parameters of G2 functions

    '''

    plt.figure(1, figsize=[7, 5])
    plt.title(r'Exemplary Visualization of $G_2$ Functions')
    plt.xlabel(r'Distance $R$ [a.u.]')
    plt.ylabel(r'$G_2(R)$')

    for rsval in rs:
        g2 = np.exp(- eta * (rr - rsval)**2) * cutoff(rr)
        plt.plot(rr, g2, label='{:.3f}'.format(rsval))

    plt.legend(title=r'$R_\mathrm{s}$ Value:')
    plt.tight_layout()
    plt.savefig('g2.svg', dpi=900, format='svg', transparent=True)
    plt.show()


def cutoff(xx):
    '''Calculates the cutoff function values for given distance array.

    Args:

        xx (1darray): contains distances to calculate the cutoff function for

    Returns:

        cut (1darray): cutoff function values

    '''

    cut = np.where(xx <= RCUT, 0.5 * (np.cos(np.pi * xx / RCUT) + 1.0), 0.0)

    return cut


if __name__ == '__main__':
    main()
