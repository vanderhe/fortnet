#!/usr/bin/env python3


'''
Application example of an H2O geometry optimization by invoking
Fortnet via the socket interface to the ASE Python package.
'''

import sys
from subprocess import Popen
from ase.io import read, write
from ase.optimize import BFGS
from ase.io.trajectory import Trajectory
from ase.calculators.socketio import SocketIOCalculator


UNIXSOCKET = 'fortnet'
FNET = 'fnet'
GEO_IN = './H2O_unoptimized.xyz'
GEO_OUT = './H2O_optimized.xyz'
TRAJECTORY = 'opt.traj'


def main():
    '''Main driver routine.'''

    system = read(GEO_IN, format='xyz')

    opt = BFGS(system, trajectory=TRAJECTORY, logfile='opt.log')

    with SocketIOCalculator(log=sys.stdout, unixsocket=UNIXSOCKET) as calc:
        Popen(FNET)
        system.set_calculator(calc)
        opt.run(fmax=1E-04)

    traj = Trajectory(TRAJECTORY)
    atoms = traj[-1]
    write(GEO_OUT, atoms, format='xyz')

if __name__ == "__main__":
    main()
