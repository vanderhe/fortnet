#!/usr/bin/env python3


'''
Application example of the Fnetout class, based on an output
file that provides the network predictions and targets, relevant
information are extracted and printed to the standard output.
'''


from fortformat import Fnetout


def main():
    '''Main driver routine.'''

    fnetout = Fnetout('fnetout.hdf5')

    mode = fnetout.mode
    print('Running mode:\n', mode)

    ndatapoints = fnetout.ndatapoints
    print('Number of datapoints in dataset:\n', ndatapoints)

    nglobaltargets = fnetout.nglobaltargets
    print('Number of system-wide targets (e.g. total energies):\n',
          nglobaltargets)

    natomictargets = fnetout.natomictargets
    print('Number of atomic targets (e.g. atomic forces):\n',
          natomictargets)

    globaltargets = fnetout.globaltargets
    print('System-wide targets:\n', globaltargets)

    atomictargets = fnetout.atomictargets
    print('Atomic targets: ', atomictargets)

    tforces = fnetout.tforces
    print('Whether atomic forces are present:\n', tforces)

    forces = fnetout.forces
    print('Atomic forces:\n', forces)

    atomicpredictions = fnetout.atomicpredictions
    print("Fortnet's predictions of atomic targets:\n", atomicpredictions)

    globalpredictions = fnetout.globalpredictions
    print("Fortnet's predictions of system-wide targets:\n", globalpredictions)

    globalpredictions_atomic = fnetout.globalpredictions_atomic
    print("Fortnet's atom-resolved predictions of system-wide targets:\n",
          globalpredictions_atomic)


if __name__ == '__main__':
    main()
