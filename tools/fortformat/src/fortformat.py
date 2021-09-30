#------------------------------------------------------------------------------#
#  FORTNET: A Behler-Parrinello-Neural-Network Implementation                  #
#  Copyright (C) 2020 - 2021  T. W. van der Heide                              #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#


'''
Basic Fortnet IO Format Classes

These basic Python classes implement the Fortnet input and output file format.
The Fnetdata class enables to create compatible HDF5 datasets, whereas the
Fnetout class extracts certain properties of the HDF5 output for later analysis.
'''


import os
import h5py
import numpy as np


# conversion factors
# (according to prog/fortnet/lib_dftbp/constants.F90)
BOHR__AA = 0.529177249
AA__BOHR = 1.0 / BOHR__AA


class Fnetdata:
    '''Basic Fortnet Input Format Class.'''


    def __init__(self, atoms=None, features=None, targets=None, atomic=False):
        '''Initializes a Fnetdata object.

        Args:

            atoms (list): list of ASE atoms objects, containing the geometries
                of the training dataset to be used
            features (list): list of numpy arrays containing external atomic
                features, used as an additional, optional input
            targets (list or 2darray): list of numpy arrays (atomic) or 2darray,
                containing the targets of the training dataset
            atomic (bool): true, if targets are atomic properties (e.g. forces)
                and false, if targets are system properties (e.g. total energy)

        '''

        self._atomic = atomic

        if atoms is None:
            self._withatoms = False
        else:
            self._withatoms = True

        if features is None:
            self._withfeatures = False
        else:
            self._withfeatures = True

        if not self._withatoms and not self._withfeatures:
            msg = 'Neither geometries nor features provided.'
            raise FnetdataError(msg)

        self._nsystems, self._ntotatoms, self._atoms, \
            self._nfeatures, self._features = \
                _checkfeatureconsistency(atoms=atoms, features=features)

        if targets is None:
            self._withtargets = False
            self._ntargets = 0
        else:
            self._withtargets = True
            self._targets, self._ntargets = \
                _checktargetconsistency(targets, self._nsystems, self._atomic)

        self._weights = np.ones((self._nsystems,), dtype=int)


    def _process_data(self):
        '''Based on the stored data, a list of dictionaries,
           containing the processed input, will be created.

        Returns:

            data (list): list of data entries for each datapoint
            unique_global_zz (1darray): unique atomic numbers of the dataset

        '''

        data = []
        tmp_zz = []

        for isys in range(self._nsystems):

            tmp = {}

            if self._atomic and self._withtargets:
                tmp['targets'] = self._targets[isys]
            elif not self._atomic and self._withtargets:
                tmp['targets'] = self._targets[isys, :]

            if self._withatoms:

                periodic = self._atoms[isys].get_pbc()

                if sum(periodic) == 3:
                    periodic = True
                elif sum(periodic) == 0:
                    periodic = False
                else:
                    msg = 'Currently, only uniform pbc are supported.'
                    raise FnetdataError(msg)

                tmp['periodic'] = periodic

                if periodic:
                    tmp['coords'] = self._atoms[isys].get_scaled_positions()
                else:
                    # dataset expects coordinates in Bohr
                    tmp['coords'] = self._atoms[isys].get_positions() * AA__BOHR

                if periodic:
                    # dataset expects lattice vectors in Bohr
                    tmp['basis'] = self._atoms[isys].get_cell()[:, :] * AA__BOHR
                else:
                    tmp['basis'] = None

                tmp['natoms'] = len(self._atoms[isys])

                tmp['atomicnumbers'] = \
                    np.array(self._atoms[isys].get_atomic_numbers())
                tmp_zz.append(set(tmp['atomicnumbers']))

                tmp['typenames'] = list(self._atoms[isys].symbols)

                # create a dictionary with unique species and id's
                atomtospecies = dict()
                localattolocalsp = np.empty(tmp['natoms'], dtype=int)
                for species in tmp['typenames']:
                    if not species in atomtospecies:
                        atomtospecies[species] = len(atomtospecies) + 1
                # mapping from local atom index to local species name
                tmp['atomtospecies'] = atomtospecies
                localtypes = list(['null'] * len(atomtospecies.keys()))
                for species in atomtospecies:
                    localtypes[atomtospecies[species] - 1] = species
                # list of string representations of local species
                tmp['localtypes'] = localtypes
                for iatom in range(tmp['natoms']):
                    ispecies = tmp['atomtospecies'][
                        tmp['typenames'][iatom]]
                    localattolocalsp[iatom] = ispecies
                # local atom index to local species index
                tmp['localattolocalsp'] = localattolocalsp

            data.append(tmp.copy())

        # extract global atomic numbers contained in the dataset
        if self._withatoms:
            unique_global_zz = {item for sublist in tmp_zz for item in sublist}
            unique_global_zz = np.sort(
                np.array(list(unique_global_zz), dtype=int))
        else:
            unique_global_zz = None

        # map local atom index to global species index
        if self._withatoms:
            for isys in range(self._nsystems):
                data[isys]['localattoglobalsp'] = \
                    np.empty(data[isys]['natoms'], dtype=int)
                for iatom in range(data[isys]['natoms']):
                    pos = np.where(
                        unique_global_zz == data[isys]['atomicnumbers'][iatom])
                    # Fortran indexing starts at 1
                    data[isys]['localattoglobalsp'][iatom] = pos[0] + 1

        return data, unique_global_zz


    def _create_contiguous_hdf(self, fname, data, zz):
        '''Creates a contiguous HDF5 dataset file.

        Args:

            fname (str): filename of dataset file to write
            data (list): dictionaries, containing the necessary information
            zz (1darray): unique atomic numbers of the dataset

        Returns:

            fid (file): in-memory hdf file, representing a contiguous dataset

        '''

        fid = h5py.File(fname, 'w')

        rootgrp = fid.create_group('fnetdata')

        datagrp = rootgrp.create_group('dataset')
        datagrp.attrs['ndatapoints'] = self._nsystems
        datagrp.attrs['nextfeatures'] = self._nfeatures
        datagrp.attrs['withstructures'] = int(self._withatoms)

        if self._withatoms:
            datagrp.attrs['ntotatoms'] = self._ntotatoms
            types = datagrp.create_dataset(
                'atomicnumbers', zz.shape, dtype='int')
            types[...] = zz

        traingrp = datagrp.create_group('training')
        traingrp.attrs['ntargets'] = self._ntargets
        traingrp.attrs['atomic'] = int(self._atomic)

        for isys in range(self._nsystems):

            subroot = datagrp.create_group('datapoint{}'.format(isys + 1))

            hdf_append_weight(subroot, self._weights[isys])

            if self._withatoms:
                hdf_append_geometry(subroot, data[isys], True)

            if self._withtargets:
                hdf_append_targets(subroot, data[isys]['targets'])

            if self._withfeatures:
                hdf_append_external_features(subroot, self._features[isys])

        return fid


    def dump(self, fname):
        '''Based on the stored data, a contiguous
           dataset file will get dumped to disk.

        Args:

            fname (str): filename of dataset file to write

        '''

        if not isinstance(fname, str):
            msg = 'Invalid dataset filename, string expected.'
            raise FnetdataError(msg)

        data, unique_global_zz = self._process_data()

        fid = self._create_contiguous_hdf(fname, data, unique_global_zz)
        dump_as_hdf(fid, fname)


    @property
    def weights(self):
        '''Defines property, providing the weight of each datapoint.

        Returns:

            weights (1darray): integer-valued array of datapoint weights

        '''

        return self._weights


    @weights.setter
    def weights(self, weights):
        '''Sets user-specified weighting of each datapoint.'''

        weights = np.array(weights)

        if weights.ndim != 1:
            msg = 'Invalid weights found, 1-dimensional list or array expected.'
            raise FnetdataError(msg)

        if not issubclass(weights.dtype.type, np.integer) or any(weights < 1):
            msg = 'Invalid weight(s) found, choose positive integers.'
            raise FnetdataError(msg)

        self._weights = weights


    @property
    def ndatapoints(self):
        '''Defines property, providing the number of datapoints.

        Returns:

            nsystems (1darray): total number of datapoints

        '''

        return self._nsystems


    @property
    def ntargets(self):
        '''Defines property, providing the number of targets.

        Returns:

            ntargets (int): if targets are atomic, the number of targets
                per atom gets returned, otherwise number of targets per system

        '''

        return self._ntargets


    @property
    def ntotatoms(self):
        '''Defines property, providing the total number of atoms in the dataset.

        Returns:

            ntotatoms (1darray): total number of atoms in the dataset

        '''

        return self._ntotatoms


def _checkfeatureconsistency(atoms=None, features=None):
    '''Performs basic consistency checks on the atomic features.

    Args:

        atoms (list): list of ASE atoms objects, containing the geometries
            of the training dataset to be used if provided
        features (list): list of numpy arrays containing external atomic
            features, used as an additional, optional input

    Returns:

        nsystems (int): number of datapoints in dataset
        ntotatoms (int): total number of atoms in geometry list
        atoms (list): list of ASE atoms objects, containing the geometries
            of the training dataset to be used if provided
        nfeatures (int): number of features per atom
        features (list): list of numpy arrays containing external atomic
            features, used as an additional, optional input

    '''

    if features is not None and not isinstance(features, list):
        msg = 'Expected external features as list.'
        raise FnetdataError(msg)

    if atoms is not None and not isinstance(atoms, list):
        msg = 'Expected geometry features as list.'
        raise FnetdataError(msg)

    if features is not None and len(features) == 0:
        msg = 'Empty list of external features provided.'
        raise FnetdataError(msg)

    if atoms is not None and len(atoms) == 0:
        msg = 'Empty list of geometry features provided.'
        raise FnetdataError(msg)

    # either geometries or external features must be present
    if atoms is not None:
        nsystems = len(atoms)
    else:
        nsystems = len(features)

    ntotatoms = 0

    if atoms is not None:
        for isys in range(nsystems):
            ntotatoms += len(atoms[isys])

    if atoms is not None and features is not None:
        msg = 'Mismatch in number of external features and ' +\
            'number of atoms of the corresponding geometry.'
        for isys in range(nsystems):
            ntotatoms += len(atoms[isys])
            if not len(atoms[isys]) == features[isys].shape[0]:
                raise FnetdataError(msg)

        nfeatures = features[0].shape[1]
    else:
        nfeatures = 0

    return nsystems, ntotatoms, atoms, nfeatures, features


def _checktargetconsistency(targets, nsystems, atomic):
    '''Performs basic consistency checks on the target values.

    Args:

        targets (list or 2darray): list of numpy arrays (atomic) or 2darray,
            containing the targets of the training dataset
        nsystems (int): number of datapoints in dataset
        atomic (bool): true, if targets are atomic properties (e.g. forces)
            and false, if targets are system properties (e.g. total energy)

    Returns:

        targets (list or 2darray): list of numpy arrays (atomic) or 2darray,
            containing the targets of the training dataset
        ntargets (int): number of features per atom

    '''

    if not atomic and isinstance(targets, list):
        targets = np.array(targets)

    if atomic and len(targets) == 0:
        msg = 'Empty list of targets provided.'
        raise FnetdataError(msg)

    if not atomic and sum(targets.shape) < 2:
        msg = 'Empty list of targets provided.'
        raise FnetdataError(msg)

    if not atomic and targets.ndim != 2:
        msg = 'Invalid number of target dimensions, ' + \
            'specify (nDatapoints, nTargets).'
        raise FnetdataError(msg)

    if atomic and nsystems != len(targets) or \
    not atomic and nsystems != targets.shape[0]:
        msg = 'Number of features and targets does not match.'
        raise FnetdataError(msg)

    if atomic:
        ntargets = targets[0].shape[1]
    else:
        ntargets = targets.shape[1]

    return targets, ntargets


def dump_as_hdf(fid, fname):
    '''Dumps a given in-memory hdf file to disk.

    Args:

        fid (file): hdf file to write to disk
        fname (str): path to write the file to

    '''

    fname = os.path.abspath(fname)
    os.makedirs(os.path.dirname(fname), exist_ok=True)

    fid.close()


def hdf_append_weight(root, weight):
    '''Appends the datapoint weight to a given in-memory hdf file.

    Args:

        root (hdf group): hdf group
        weight (int): positive integer weight of current datapoint

    '''

    root.attrs['weight'] = weight


def hdf_append_geometry(root, data, frac):
    '''Appends geometry information to a given in-memory hdf file.

    Args:

        root (hdf group): hdf group
        data (dict): dictionary, containing the necessary information
        frac (bool): true, if coordinates should be stored in units of the
            lattice vectors (presupposes a periodic structure)

    '''

    geogrp = root.create_group('geometry')

    geogrp.attrs['fractional'] = int(frac and data['periodic'])
    geogrp.attrs['localtypes'] = ','.join(data['localtypes'])

    localattolocalsp = geogrp.create_dataset(
        'localattolocalsp', data['localattolocalsp'].shape, dtype='int')
    localattolocalsp[...] = data['localattolocalsp']

    localattoglobalsp = geogrp.create_dataset(
        'localattoglobalsp', data['localattoglobalsp'].shape, dtype='int')
    localattoglobalsp[...] = data['localattoglobalsp']

    localattoatnum = geogrp.create_dataset(
        'localattoatnum', data['atomicnumbers'].shape, dtype='int')
    localattoatnum[...] = data['atomicnumbers']

    coords = geogrp.create_dataset(
        'coordinates', data['coords'].shape, dtype='float')
    coords[...] = data['coords']

    geogrp.attrs['periodic'] = int(data['periodic'])

    if data['periodic']:
        basis = geogrp.create_dataset(
            'basis', data['basis'].shape, dtype='float')
        basis[...] = data['basis']


def hdf_append_external_features(root, data):
    '''Appends external atomic features to a given in-memory hdf file.

    Args:

        root (hdf group): hdf group
        data (2darray): numpy array containing external atomic
            features, used as an additional, optional input

    '''

    features = root.create_dataset('extfeatures', data.shape, dtype='float')
    features[...] = data


def hdf_append_targets(root, data):
    '''Appends target information to a given in-memory hdf file.

    Args:

        root (hdf group): hdf group
        data (2darray): atomic or global targets

    '''

    if data.ndim == 1:
        tmp = np.empty((len(data), 1), dtype=float)
        tmp[:, 0] = data
    else:
        tmp = data

    targets = root.create_dataset('targets', tmp.shape, dtype='float')
    targets[...] = tmp


class Fnetout:
    '''Basic Fortnet Output Format Class.'''


    def __init__(self, fname):
        '''Initializes a Fnetout object.

        Args:

            fname (str): filename to extract data from

        '''

        self._fname = fname

        with h5py.File(self._fname, 'r') as fnetoutfile:
            fnetout = fnetoutfile['fnetout']
            self._mode = fnetout.attrs.get('mode').decode('UTF-8').strip()
            if not self._mode in ('validate', 'predict'):
                raise FnetoutError('Invalid running mode specification.')

            output = fnetoutfile['fnetout']['output']
            self._ndatapoints = output.attrs.get('ndatapoints')
            if len(self._ndatapoints) == 1:
                # number of datapoints stored in array of size 1
                self._ndatapoints = self._ndatapoints[0]
            else:
                msg = "Error while reading fnetout file '" + self._fname + \
                    "'. Unrecognized number of datapoints obtained."
                raise FnetoutError(msg)
            self._tforces = output.attrs.get('tforces')
            if len(self._tforces) == 1:
                # booleans stored in integer arrays of size 1
                self._tforces = bool(self._tforces[0])
            else:
                msg = "Error while reading fnetout file '" + self._fname + \
                    "'. Unrecognized force specification obtained."
                raise FnetoutError(msg)

            self._targettype = \
                output.attrs.get('targettype').decode('UTF-8').strip()
            if not self._targettype in ('atomic', 'global'):
                raise FnetoutError('Invalid running mode obtained.')

            # get number of atomic or global predictions/targets
            self._npredictions = np.shape(
                np.array(output['datapoint1']['output']))[1]

            if self._mode == 'validate':
                self._npredictions = int(self._npredictions / 2)


    @property
    def mode(self):
        '''Defines property, providing the mode of the Fortnet run.

        Returns:

            mode (str): mode of the run that produced the Fnetout file

        '''

        return self._mode


    @property
    def ndatapoints(self):
        '''Defines property, providing the number of datapoints.

        Returns:

            ndatapoints (int): total number of datapoints of the training

        '''

        return self._ndatapoints


    @property
    def targettype(self):
        '''Defines property, providing the target type.

        Returns:

            targettype (str): type of targets the network was trained on

        '''

        return self._targettype


    @property
    def tforces(self):
        '''Defines property, providing hint whether atomic forces are present.

        Returns:

            tforces (bool): true, if atomic forces are supplied

        '''

        return self._tforces


    @property
    def predictions(self):
        '''Defines property, providing the predictions of Fortnet.

        Returns:

            predictions (list or 2darray): predictions of the network

        '''

        with h5py.File(self._fname, 'r') as fnetoutfile:
            output = fnetoutfile['fnetout']['output']
            if self._targettype == 'atomic':
                predictions = []
                for idata in range(self._ndatapoints):
                    dataname = 'datapoint' + str(idata + 1)
                    if self._mode == 'validate':
                        predictions.append(
                            np.array(output[dataname]['output'],
                                     dtype=float)[:, :self._npredictions])
                    else:
                        predictions.append(
                            np.array(output[dataname]['output'], dtype=float))
            else:
                predictions = np.empty(
                    (self._ndatapoints, self._npredictions), dtype=float)
                for idata in range(self._ndatapoints):
                    dataname = 'datapoint' + str(idata + 1)
                    if self._mode == 'validate':
                        predictions[idata, :] = \
                            np.array(output[dataname]['output'],
                                     dtype=float)[0, :self._npredictions]
                    else:
                        predictions[idata, :] = \
                            np.array(output[dataname]['output'],
                                     dtype=float)[0, :]

        return predictions


    @property
    def targets(self):
        '''Defines property, providing the targets during training.

        Returns:

            targets (list or 2darray): targets during training

        '''

        if self._mode == 'predict':
            return None

        with h5py.File(self._fname, 'r') as fnetoutfile:
            output = fnetoutfile['fnetout']['output']
            if self._targettype == 'atomic':
                targets = []
                for idata in range(self._ndatapoints):
                    dataname = 'datapoint' + str(idata + 1)
                    targets.append(
                        np.array(output[dataname]['output'],
                                 dtype=float)[:, self._npredictions:])
            else:
                targets = np.empty(
                    (self._ndatapoints, self._npredictions), dtype=float)
                for idata in range(self._ndatapoints):
                    dataname = 'datapoint' + str(idata + 1)
                    targets[idata, :] = \
                        np.array(output[dataname]['output'],
                                 dtype=float)[0, self._npredictions:]

        return targets


    @property
    def forces(self):
        '''Defines property, providing the atomic forces, if supplied.

        Returns:

            forces (list): atomic forces on atoms

        '''

        tmp1 = []

        if self._targettype == 'atomic':
            msg = "Error while extracting forces from fnetout file '" \
                + self._fname + \
                "'. Forces only supplied for global property targets."
            raise FnetoutError(msg)

        with h5py.File(self._fname, 'r') as fnetoutfile:
            output = fnetoutfile['fnetout']['output']
            for idata in range(self._ndatapoints):
                dataname = 'datapoint' + str(idata + 1)
                tmp1.append(np.array(output[dataname]['forces'], dtype=float))

        # convert to shape np.shape(forces[iData][iTarget]) = (iAtom, 3)
        forces = []
        for tmp2 in tmp1:
            entry = []
            if not np.shape(tmp2)[1]%3 == 0:
                msg = "Error while extracting forces from fnetout file '" \
                    + self._fname + \
                    "'. Expected three force components and global target."
                raise FnetoutError(msg)
            for jj in range(int(np.shape(tmp2)[1] / 3)):
                entry.append(tmp2[:, 3 * jj:3 * (jj + 1)])
            forces.append(entry)

        return forces


class FnetdataError(Exception):
    '''Exception thrown by the Fnetdata class.'''


class FnetoutError(Exception):
    '''Exception thrown by the Fnetout class.'''
