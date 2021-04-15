#------------------------------------------------------------------------------#
#  FORTNET: A Behler-Parrinello-Neural-Network Implementation                  #
#  Copyright (C) 2020 - 2021  T. W. van der Heide                              #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#


'''
Basic Fortnet Input Format Class

This basic Python class implements the Fortnet input file format. The input
geometries, stored in a list of ASE atoms objects, and targets, stored in a
list of Numpy arrays, conveniently get dumped to disk as simple text files.
'''


import os
import xml.etree.ElementTree as ET
import numpy as np


# conversion factors
# (according to prog/fortnet/lib_dftbp/constants.F90)
BOHR__AA = 0.529177249
AA__BOHR = 1.0 / BOHR__AA


class Fortformat:
    '''Basic Fortnet Input Format Class.'''


    def __init__(self, atoms, paths, targets=None, atomic=False, frac=False):
        '''Initializes a Fortformat object.

        Args:

            atoms (list): list of ASE atoms objects, containing the geometries
                of the training dataset to be used
            paths (list) : list of strings, containing the paths of the output
                fnetdata.xml files, which will contain all relevant informations
            targets (list or 2darray): list of numpy arrays (atomic) or 2darray,
                containing the targets of the training dataset
            atomic (bool): true, if targets are atomic properties (e.g. forces)
                and false, if targets are system properties (e.g. total energy)
            frac (bool): true, coordinates should be stored in units of the
                lattice vectors (presupposes a periodic structure)

        '''

        if targets is not None:
            self._withtargets = True
        else:
            self._withtargets = False

        self._paths = paths
        self._atomic = atomic
        self._frac = frac

        self.checkfeatureconsistency(atoms)
        self._atoms = atoms
        self._nsystems = len(self._atoms)

        if self._withtargets:
            self.checktargetconsistency(targets)
            self._targets = targets
        else:
            self._targets = None
            self._ntargets = None


    def dump(self):
        '''Based on the stored data, fnetdata.xml files get dumped to disk.'''

        data = {}

        for isys in range(self._nsystems):

            periodic = self._atoms[isys].get_pbc()

            if sum(periodic) == 3:
                periodic = True
            elif sum(periodic) == 0:
                periodic = False
            else:
                msg = 'Currently, only uniform pbc are supported.'
                raise FortformatError(msg)

            data['periodic'] = periodic

            if self._frac and periodic:
                data['coords'] = self._atoms[isys].get_scaled_positions()
            else:
                # fnetdata.xml expects coordinates in Bohr
                data['coords'] = self._atoms[isys].get_positions() * AA__BOHR

            if periodic:
                # fnetdata.xml expects lattice vectors in Bohr
                data['basis'] = self._atoms[isys].get_cell()[:, :] * AA__BOHR
            else:
                data['basis'] = None

            data['natoms'] = len(self._atoms[isys])

            if self._atomic and self._withtargets:
                data['targets'] = self._targets[isys]
            elif not self._atomic and self._withtargets:
                data['targets'] = self._targets[isys, :]

            data['typenames'] = list(self._atoms[isys].symbols)

            # create a dictionary with unique species and id's
            atomtospecies = dict()
            for species in data['typenames']:
                if not species in atomtospecies:
                    atomtospecies[species] = len(atomtospecies) + 1
            uniquespecies = list(['null'] * len(atomtospecies.keys()))
            for species in atomtospecies:
                uniquespecies[atomtospecies[species] - 1] = species
            data['atomtospecies'] = atomtospecies
            data['uniquespecies'] = uniquespecies
            data['uniquespecies'] = ['"' + entry + '"' for entry in
                                     data['uniquespecies']]

            if not self._paths[isys].endswith('.xml'):
                filename = os.path.join(self._paths[isys], 'fnetdata.xml')
            else:
                filename = self._paths[isys]

            self.dump_as_xml(data, filename)


    def dump_as_xml(self, data, filename):
        '''Dumps a single fnetdata.xml file to disk.

        Args:
            data (dict): dictionary, containing all the necessary informations
            filename (str): path where to write output file to

        '''

        nl = '\n'
        floatfmt = '{:24.16e}'

        root = ET.Element('fnetdata')

        geotag = ET.SubElement(root, 'geometry')

        tmptag = ET.SubElement(geotag, 'typenames')
        tmptag.text = ' '.join(data['uniquespecies'])

        tmptag = ET.SubElement(geotag, 'fractional')
        if self._frac:
            tmptag.text = 'Yes'
        else:
            tmptag.text = 'No'

        tmptag = ET.SubElement(geotag, 'typesandcoordinates')

        tmp = []
        for iatom in range(data['natoms']):
            ispecies = data['atomtospecies'][data['typenames'][iatom]]
            coordentry = data['coords'][iatom, :]
            tmp.append('  {} '.format(ispecies) + \
                       (3 * floatfmt).format(*coordentry))
        tmptag.text = nl + '\n'.join(tmp) + nl

        tmptag = ET.SubElement(geotag, 'periodic')
        if data['periodic']:
            tmptag.text = 'Yes'
        else:
            tmptag.text = 'No'

        tmptag = ET.SubElement(geotag, 'latticevectors')

        tmp = []
        for ilatt in range(3):
            tmp.append((3 * floatfmt).format((*data['basis'][ilatt, :])))
        tmptag.text = nl + '\n'.join(tmp) + nl

        tmptag = ET.SubElement(geotag, 'coordinateorigin')
        tmptag.text = nl + \
            (3 * floatfmt).format(*np.array((0.0, 0.0, 0.0))) + nl

        if self._withtargets:

            traintag = ET.SubElement(root, 'training')

            tmptag = ET.SubElement(traintag, 'atomic')
            if self._atomic:
                tmptag.text = 'Yes'
            else:
                tmptag.text = 'No'

            tmptag = ET.SubElement(traintag, 'ntargets')
            tmptag.text = '{}'.format(self._ntargets)

            tmptag = ET.SubElement(traintag, 'targets')
            tmp = []
            if self._atomic:
                for iatom in range(data['natoms']):
                    tmp.append((self._ntargets * \
                                   floatfmt).format(*data['targets'][iatom, :]))
                tmptag.text = nl + '\n'.join(tmp) + nl
            else:
                tmptag.text = nl + (self._ntargets * \
                    floatfmt).format(*data['targets']) + nl

        fnetdata = ET.ElementTree(_beautify(root))
        fnetdata.write(filename, xml_declaration=True, encoding='utf-8')


    def checktargetconsistency(self, targets):
        '''Performs basic consistency checks on the target values.

        Args:

            targets (list or 2darray): list of numpy arrays (atomic) or 2darray,
                containing the targets of the training dataset

        '''

        if not self._atomic and isinstance(targets, list):
            targets = np.array(targets)

        if self._atomic and len(targets) == 0:
            msg = 'Empty list of targets provided.'
            raise FortformatError(msg)

        if not self._atomic and sum(targets.shape) < 2:
            msg = 'Empty list of targets provided.'
            raise FortformatError(msg)

        if not self._atomic and targets.ndim != 2:
            msg = 'Invalid number of target dimensions, ' + \
                'specify (nDatapoints, nTargets).'
            raise FortformatError(msg)

        if self._atomic and self._nsystems != len(targets) or \
        not self._atomic and self._nsystems != targets.shape[0]:
            msg = 'Number of features and targets does not match.'
            raise FortformatError(msg)

        if self._atomic:
            self._ntargets = targets[0].shape[1]
        else:
            self._ntargets = targets.shape[1]


    def checkfeatureconsistency(self, atoms):
        '''Performs basic consistency checks on the atomic features.

        Args:

            atoms (list): list of ASE atoms objects, containing the geometries
                of the training dataset to be used

        '''

        if not isinstance(atoms, list):
            msg = 'Expected input features as list.'
            raise FortformatError(msg)

        if len(atoms) == 0:
            msg = 'Empty list of features provided.'
            raise FortformatError(msg)

        if len(self._paths) != len(atoms):
            msg = 'Number of features and paths does not match.'
            raise FortformatError(msg)


    def get_nr_structures(self):
        '''Queries the number of structures.

        Returns:

            self._nsystems (int): number of structures

        '''

        return self._nsystems


    def get_nr_targets(self):
        '''Queries the number of targets.

        Returns:

            self._ntargets (int): If targets are atomic, the number of targets
                per atom gets returned (otherwise number of targets per system).

        '''

        return self._ntargets


def _beautify(root, level=0, nindent=0):
    '''Break lines and indent xml tags, if desired.'''

    indent = nindent * ' '

    ii = '\n' + level * indent
    jj = '\n' + (level - 1) * indent

    if len(root) > 0:
        if not root.text or not root.text.strip():
            root.text = ii + indent
        if not root.tail or not root.tail.strip():
            root.tail = ii
        for subroot in root:
            _beautify(subroot, level + 1)
        if not root.tail or not root.tail.strip():
            root.tail = jj
    else:
        if level and (not root.tail or not root.tail.strip()):
            root.tail = jj

    return root


class FortformatError(Exception):
    '''Exception thrown by the Fortformat class.'''
