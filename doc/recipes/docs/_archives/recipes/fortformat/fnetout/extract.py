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
    print('Running mode: ', mode)

    ndatapoints = fnetout.ndatapoints
    print('Number of datapoints in training: ', ndatapoints)

    targettype = fnetout.targettype
    print('Type of targets: ', targettype)

    predictions = fnetout.predictions
    print("Fortnet's predictions: ", predictions)

    targets = fnetout.targets
    print('Targets while trained: ', targets)


if __name__ == '__main__':
    main()
