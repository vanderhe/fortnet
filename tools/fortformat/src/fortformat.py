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

NL = '\n'
FLOATFMT = '{:24.16e}'


class Fortformat:
    '''Basic Fortnet Input Format Class.'''


    def __init__(self, atoms, paths, targets=None, atomic=False,
                 frac=False):
        '''Initializes a Fortformat object.

        Args:

            atoms (list): list of ASE atoms objects, containing the geometries
                of the training dataset to be used
            paths (list) : string or list of strings, containing a single path
                to a contiguous dataset file or multiple paths to single
                fnetdata.xml files, which will contain all relevant information
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

        if isinstance(paths, list):
            self._contiguous = False
        elif isinstance(paths, str):
            self._contiguous = True
        else:
            msg = \
                'Invalid format of paths provided. Choose either a single ' +\
                'path for a contiguous dataset file or a list of paths for ' +\
                'single fnetdata.xml files.'
            raise FortformatError(msg)

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

        self.process_paths()


    def process_paths(self):
        '''In case of non-contiguous files, appends suffix if not present.'''

        if not self._contiguous:
            for isys in range(self._nsystems):
                if not self._paths[isys].endswith('.xml'):
                    self._paths[isys] = os.path.join(
                        self._paths[isys], 'fnetdata.xml')


    def process_data(self):
        '''Based on the stored data, a list of dictionaries,
           containing the processed input, will be created.
        '''

        data = []
        tmp = {}

        for isys in range(self._nsystems):

            periodic = self._atoms[isys].get_pbc()

            if sum(periodic) == 3:
                periodic = True
            elif sum(periodic) == 0:
                periodic = False
            else:
                msg = 'Currently, only uniform pbc are supported.'
                raise FortformatError(msg)

            tmp['periodic'] = periodic

            if self._frac and periodic:
                tmp['coords'] = self._atoms[isys].get_scaled_positions()
            else:
                # fnetdata.xml expects coordinates in Bohr
                tmp['coords'] = self._atoms[isys].get_positions() * AA__BOHR

            if periodic:
                # fnetdata.xml expects lattice vectors in Bohr
                tmp['basis'] = self._atoms[isys].get_cell()[:, :] * AA__BOHR
            else:
                tmp['basis'] = None

            tmp['natoms'] = len(self._atoms[isys])

            if self._atomic and self._withtargets:
                tmp['targets'] = self._targets[isys]
            elif not self._atomic and self._withtargets:
                tmp['targets'] = self._targets[isys, :]

            tmp['typenames'] = list(self._atoms[isys].symbols)

            # create a dictionary with unique species and id's
            atomtospecies = dict()
            for species in tmp['typenames']:
                if not species in atomtospecies:
                    atomtospecies[species] = len(atomtospecies) + 1
            uniquespecies = list(['null'] * len(atomtospecies.keys()))
            for species in atomtospecies:
                uniquespecies[atomtospecies[species] - 1] = species
            tmp['atomtospecies'] = atomtospecies
            tmp['uniquespecies'] = uniquespecies
            tmp['uniquespecies'] = ['"' + entry + '"' for entry in
                                     tmp['uniquespecies']]

            data.append(tmp.copy())

        return data


    def create_single_xml(self, data):
        '''Creates a single fnetdata.xml file, as ET xml-tree instance.

        Args:

            data (dict): dictionary, containing all the necessary information

        Returns:

            root (ET instance): xml-tree, representing a single dataset file

        '''

        root = ET.Element('fnetdata')

        root = xml_append_geometry(root, data, self._frac)

        if self._withtargets:
            root = xml_append_targets(root, data, self._ntargets,
                                      self._atomic, contiguous=False)

        return root


    def create_contiguous_xml(self, data):
        '''Creates a contiguous fnetdata.xml file, as ET xml-tree instance.

        Args:

            data (list): dictionaries, containing the necessary information

        Returns:

            root (ET instance): xml-tree, representing a contiguous dataset file

        '''

        root = ET.Element('fnetdata')

        datasettag = ET.SubElement(root, 'dataset')
        tmptag = ET.SubElement(datasettag, 'ndatapoints')
        tmptag.text = '{}'.format(self._nsystems)

        if self._withtargets:

            traintag = ET.SubElement(root, 'training')

            tmptag = ET.SubElement(traintag, 'atomic')
            if self._atomic:
                tmptag.text = 'Yes'
            else:
                tmptag.text = 'No'

            tmptag = ET.SubElement(traintag, 'ntargets')
            tmptag.text = '{}'.format(self._ntargets)

        for isys in range(self._nsystems):

            subroot = ET.SubElement(root, 'datapoint{}'.format(isys + 1))

            subroot = xml_append_geometry(subroot, data[isys], self._frac)

            if self._withtargets:
                subroot = xml_append_targets(
                    subroot, data[isys], self._ntargets,
                    self._atomic, contiguous=True)

        return root


    def dump(self):
        '''Based on the stored data, either a contiguous fnetdata
           or various single fnetdata files get dumped to disk.
        '''

        data = self.process_data()

        if self._contiguous:
            xml = self.create_contiguous_xml(data)
            dump_as_xml(xml, self._paths)
        else:
            for isys in range(self._nsystems):
                xml = self.create_single_xml(data[isys])
                dump_as_xml(xml, self._paths[isys])


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

        if not self._contiguous and len(self._paths) != len(atoms):
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


def dump_as_xml(root, fname):
    '''Dumps a given xml-tree to disk.

    Args:

        root (ET instance): xml-tree to write to disk
        fname (str): path to write the file to

    '''

    fname = os.path.abspath(fname)
    os.makedirs(os.path.dirname(fname), exist_ok=True)

    fnetdata = ET.ElementTree(_beautify(root))
    fnetdata.write(fname, xml_declaration=True, encoding='utf-8')


def xml_append_geometry(root, data, frac):
    '''Appends geometry information to a given xml-tree root.

    Args:

        root (ET instance): xml-tree
        data (dict): dictionary, containing the necessary information
        frac (bool): true, if coordinates should be stored in units of the
            lattice vectors (presupposes a periodic structure)

    Returns:

        root (ET instance): xml-tree with a geometry appended to it

    '''

    geotag = ET.SubElement(root, 'geometry')

    tmptag = ET.SubElement(geotag, 'typenames')
    tmptag.text = ' '.join(data['uniquespecies'])

    tmptag = ET.SubElement(geotag, 'fractional')
    if frac:
        tmptag.text = 'Yes'
    else:
        tmptag.text = 'No'

    tmptag = ET.SubElement(geotag, 'typesandcoordinates')

    tmp = []
    for iatom in range(data['natoms']):
        ispecies = data['atomtospecies'][
            data['typenames'][iatom]]
        coordentry = data['coords'][iatom, :]
        tmp.append('  {} '.format(ispecies) + \
                   (3 * FLOATFMT).format(*coordentry))
    tmptag.text = NL + '\n'.join(tmp) + NL

    tmptag = ET.SubElement(geotag, 'periodic')
    if data['periodic']:
        tmptag.text = 'Yes'
    else:
        tmptag.text = 'No'

    tmptag = ET.SubElement(geotag, 'latticevectors')

    tmp = []
    for ilatt in range(3):
        tmp.append((3 * FLOATFMT).format((*data['basis'][ilatt, :])))
    tmptag.text = NL + '\n'.join(tmp) + NL

    tmptag = ET.SubElement(geotag, 'coordinateorigin')
    tmptag.text = NL + \
        (3 * FLOATFMT).format(*np.array((0.0, 0.0, 0.0))) + NL

    return root


def xml_append_targets(root, data, ntargets, atomic, contiguous=False):
    '''Appends target information to a given xml-tree root.

    Args:

        root (ET instance): xml-tree
        data (dict): dictionary, containing the necessary information
        ntargets (int): number of atomic targets if atomic is true,
            otherwise number of global system targets
        atomic (bool): true, if targets are atomic properties (e.g. forces)
            and false, if targets are system properties (e.g. total energy)
        contiguous (bool): decides whether reduced target information
            (contiguous format) or complete entries (single format) should get
            appended (default: False)

    Returns:

        root (ET instance): xml-tree with targets appended to it

    '''

    tmp = []

    traintag = ET.SubElement(root, 'training')

    if not contiguous:
        tmptag = ET.SubElement(traintag, 'atomic')
        if atomic:
            tmptag.text = 'Yes'
        else:
            tmptag.text = 'No'

        tmptag = ET.SubElement(traintag, 'ntargets')
        tmptag.text = '{}'.format(ntargets)

    tmptag = ET.SubElement(traintag, 'targets')

    if atomic:
        for iatom in range(data['natoms']):
            tmp.append((ntargets * FLOATFMT)
                       .format(*data['targets'][iatom, :]))
        tmptag.text = NL + '\n'.join(tmp) + NL
    else:
        tmptag.text = NL + (ntargets * \
            FLOATFMT).format(*data['targets']) + NL

    return root


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
