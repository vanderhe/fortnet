#------------------------------------------------------------------------------#
#  FORTNET: A Behler-Parrinello-Neural-Network Implementation                  #
#  Copyright (C) 2020 - 2021  T. W. van der Heide                              #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#


'''
Basic Fortnet Regression Test Class

Initialized with the path to a working directory this class figures out the
rudimentary specifications of the Fortnet run and compares the existing
reference files with the output of the current version of the program.
'''


import os
import warnings
import xml.etree.ElementTree as ET
import numpy as np


ATOL=1e-06
RTOL=1e-05

ACSFOUT = 'acsf.out'
PRECOUT = 'precond.out'
FNETOUT = 'fnetout.xml'


class TestWithWorkDir:
    '''Basic Fortnet Regression Test Class.'''


    def __init__(self, wdir):
        '''Initializes a TestWithWorkDir object.

        Args:

            wdir (str): path to working directory to perform regression testing

        '''

        self._wdir = wdir

        # check for nessecary output files
        if not os.path.isfile(os.path.join(self._wdir, '_' + ACSFOUT)):
            raise TestWithWorkDirError('Reference acsf output file absent.')

        # try to figure out the expected configuration of the Fortnet run
        self._isfnetout = os.path.isfile(
            os.path.join(self._wdir, '_' + FNETOUT))
        self._isprec = os.path.isfile(os.path.join(self._wdir, '_' + PRECOUT))

        if self._isfnetout and \
           not os.path.isfile(os.path.join(self._wdir, FNETOUT)):
            raise TestWithWorkDirError('Fnetout xml file absent.')

        if self._isprec and \
           not os.path.isfile(os.path.join(self._wdir, PRECOUT)):
            raise TestWithWorkDirError('Preconditioning file absent.')

        if not self._isprec and \
           os.path.isfile(os.path.join(self._wdir, PRECOUT)):
            raise TestWithWorkDirError('Unexpected preconditioning file found.')

        # get all netstat file names
        self._refnetfiles = []
        self._refnetfiles += [fname for fname in os.listdir(self._wdir) if
                              fname.startswith('_') and fname.endswith('.net')]

        if len(self._refnetfiles) == 0:
            raise TestWithWorkDirError('No netstat files found.')

        # generate corresponding filenames to expect for testing
        self._curnetfiles = [entry.strip('_') for entry in self._refnetfiles]


    def test(self):
        '''Performs regression testing by comparing all relevant files.

        Returns:

            passed (bool): True, if all tests have passed, otherwise False

        '''

        curacsf = Acsf(os.path.join(self._wdir, ACSFOUT))
        refacsf = Acsf(os.path.join(self._wdir, '_' + ACSFOUT))
        acsfpassed = curacsf.equals(refacsf)

        if self._isfnetout:
            curfnetout = Fnetout(os.path.join(self._wdir, FNETOUT))
            reffnetout = Fnetout(os.path.join(self._wdir, '_' + FNETOUT))
            fnetoutpassed = curfnetout.equals(reffnetout)
        else:
            fnetoutpassed = True

        if self._isprec:
            curprec = Precond(os.path.join(self._wdir, PRECOUT))
            refprec = Precond(os.path.join(self._wdir, '_' + PRECOUT))
            precpassed = curprec.equals(refprec)
        else:
            precpassed = True

        singlenetpassed = []

        for curnetfile, refnetfile in zip(self._curnetfiles, self._refnetfiles):
            refnet = Netfile(os.path.join(self._wdir, refnetfile))
            curnet = Netfile(os.path.join(self._wdir, curnetfile))
            singlenetpassed.append(curnet.equals(refnet))

        netspassed = all(singlenetpassed)

        passed = acsfpassed and netspassed and precpassed and fnetoutpassed

        return passed


class Netfile:
    '''Fortnet network status output instance.'''


    def __init__(self, fname):
        '''Initializes a Netfile object.

        Args:

            fname (str): path to working directory to perform regression testing

        '''

        self._fname = fname

        with open(self._fname, 'r') as infile:
            content = infile.readlines()

        self.nettype = content[0].strip().split()[0].lower()
        self.targettype = content[0].strip().split()[1].lower()
        self.species = content[1].strip().lower()

        if int(content[2].strip()) >= 2:
            self.nlayer = int(content[2].strip())
        else:
            msg = 'Total number of layers must be greater or equal to two.'
            raise NetfileError(msg)

        self.dims = np.array(content[3].strip().split(), dtype=int)

        if any([dim <= 0 for dim in self.dims]):
            msg = 'There must be at least one neuron per layer.'
            raise NetfileError(msg)

        self.activation = content[4].strip().lower()

        nbiases = np.sum(self.dims) - self.dims[0]
        self.biases = np.loadtxt(self._fname, skiprows=5, max_rows=nbiases)

        self.weights = np.loadtxt(self._fname, skiprows=5+nbiases)


    def equals(self, ref, atol=ATOL, rtol=RTOL):
        '''Checks equality with another reference instance

        Args:

            ref (Acsf): reference instance to compare with
            atol (float): required absolute tolerance
            rtol (float): required relative tolerance

        Returns:

            equal (bool): True, if the two instances are equal within tolerance

        '''

        equal = self.nettype == ref.nettype

        if not equal:
            warnings.warn('Mismatch in network type.')
            return False

        equal = self.targettype == ref.targettype

        if not equal:
            warnings.warn('Mismatch in target type.')
            return False

        equal = self.species == ref.species

        if not equal:
            warnings.warn('Mismatch in species type.')
            return False

        equal = self.activation == ref.activation

        if not equal:
            warnings.warn('Mismatch in transfer function specification.')
            return False

        equal = self.nlayer == ref.nlayer

        if not equal:
            warnings.warn('Mismatch in total number of layers.')
            return False

        equal = np.allclose(self.dims, ref.dims, rtol=rtol, atol=atol)

        if not equal:
            warnings.warn('Mismatch in number of neurons per layer.')
            return False

        equal = np.allclose(self.biases, ref.biases, rtol=rtol, atol=atol)

        if not equal:
            warnings.warn('Mismatch in bias parameters.')
            return False

        equal = np.allclose(self.weights, ref.weights, rtol=rtol, atol=atol)

        if not equal:
            warnings.warn('Mismatch in weight parameters.')
            return False

        return True


class Acsf:
    '''Atom-Centered-Symmetry-Function output instance.'''


    def __init__(self, fname):
        '''Initializes a Acsf object.

        Args:

            fname (str): path to working directory to perform regression testing

        '''

        self._fname = fname

        with open(self._fname, 'r') as infile:
            content = infile.readlines()

        if float(content[0]) > 0.0:
            self.rcut = float(content[0])
        else:
            msg = 'Cutoff radius must be positive.'
            raise AcsfError(msg)

        if int(content[1].strip().split()[0]) >= 0:
            self.nradial = int(content[1].strip().split()[0])
            self.g2eta = float(content[1].strip().split()[1])
        else:
            msg = 'Number of radial functions must be greater or equal zero.'
            raise AcsfError(msg)

        if int(content[self.nradial + 2].strip().split()[0]) >= 0:
            self.nangular = int(content[self.nradial + 2].strip().split()[0])
            self.g4eta = float(content[self.nradial + 2].strip().split()[1])
        else:
            msg = 'Number of angular functions must be greater or equal zero.'
            raise AcsfError(msg)

        self.g2rs = np.empty(self.nradial)
        self.g4xilambda = np.empty((self.nangular, 2))

        for ii in range(self.nradial):
            self.g2rs[ii] = float(content[ii + 2].strip().split()[0])

        for ii in range(self.nangular):
            self.g4xilambda[ii, 0] = float(
                content[ii + self.nradial + 3].strip().split()[0])
            self.g4xilambda[ii, 1] = float(
                content[ii + self.nradial + 3].strip().split()[1])

        if int(content[self.nradial + self.nangular + 3].strip().split()[0]) \
           > 0:
            self.nspecies = int(content[self.nradial + self.nangular + 3]
                                .strip().split()[0])
        else:
            msg = 'Number of atomic species must be greater zero.'
            raise AcsfError(msg)

        self.speciesnames = []
        self.speciesids = np.empty(self.nspecies)
        for ii in range(self.nspecies):
            self.speciesnames.append(str(content[
                ii + self.nradial + self.nangular + 4].strip().split()[0]))
            self.speciesids[ii] = float(content[
                ii + self.nradial + self.nangular + 4].strip().split()[1])

        self.isprec = int(
            content[self.nradial + self.nangular + self.nspecies + 4]
            .strip().split()[0])

        if bool(self.isprec):
            self.prec = np.empty((self.nradial + self.nangular, 2))
            for ii in range(self.nradial + self.nangular):
                self.prec[ii, 0] = float(content[
                    ii + self.nradial + self.nangular + self.nspecies + 5]
                                         .strip().split()[0])
                self.prec[ii, 1] = float(content[
                    ii + self.nradial + self.nangular + self.nspecies + 5]
                                         .strip().split()[1])


    def equals(self, ref, atol=ATOL, rtol=RTOL):
        '''Checks equality with another reference instance

        Args:

            ref (Acsf): reference instance to compare with
            atol (float): required absolute tolerance
            rtol (float): required relative tolerance

        Returns:

            equal (bool): True, if the two instances are equal within tolerance

        '''

        equal = np.allclose([self.rcut], [ref.rcut], rtol=rtol, atol=atol)

        if not equal:
            warnings.warn('Mismatch in cutoff radii.')
            return False

        equal = self.nradial == ref.nradial

        if not equal:
            warnings.warn('Mismatch in number of radial mappings.')
            return False

        equal = self.nangular == ref.nangular

        if not equal:
            warnings.warn('Mismatch in number of angular mappings.')
            return False

        equal = np.allclose(self.g2rs, ref.g2rs, rtol=rtol, atol=atol)

        if not equal:
            warnings.warn('Mismatch in gaussian shift parameters.')
            return False

        equal = np.allclose([self.g2eta], [ref.g2eta], rtol=rtol, atol=atol)

        if not equal:
            warnings.warn('Mismatch in radial eta parameter.')
            return False

        equal = np.allclose([self.g4eta], [ref.g4eta], rtol=rtol, atol=atol)

        if not equal:
            warnings.warn('Mismatch in angular eta parameter.')
            return False

        equal = np.allclose(self.g4xilambda, ref.g4xilambda,
                            rtol=rtol, atol=atol)

        if not equal:
            warnings.warn('Mismatch in angular xi or lambda parameters.')
            return False

        equal = self.isprec == ref.isprec

        if not equal:
            warnings.warn('Mismatch in preconditioning specification.')
            return False

        if self.isprec and equal:
            equal = np.allclose(self.prec, ref.prec, rtol=rtol, atol=atol)
            if not equal:
                warnings.warn('Mismatch in preconditioning parameters.')
                return False

        return True


class Precond:
    '''Preconditioning output instance.'''


    def __init__(self, fname):
        '''Initializes a Precond object.

        Args:

            fname (str): path to working directory to perform regression testing

        '''

        self._fname = fname

        with open(self._fname, 'r') as infile:
            content = infile.readlines()

        if int(content[0].strip()) > 0:
            self.ntargets = int(content[0].strip())
        else:
            msg = 'Number of preconditioning tuples must be positive.'
            raise PrecondError(msg)

        self.prec = np.empty((self.ntargets, 2))

        for ii in range(self.ntargets):
            self.prec[ii, 0] = float(content[ii + 1].strip().split()[0])
            self.prec[ii, 1] = float(content[ii + 1].strip().split()[1])


    def equals(self, ref, atol=ATOL, rtol=RTOL):
        '''Checks equality with another reference instance

        Args:

            ref (Acsf): reference instance to compare with
            atol (float): required absolute tolerance
            rtol (float): required relative tolerance

        Returns:

            equal (bool): True, if the two instances are equal within tolerance

        '''

        equal = self.ntargets == ref.ntargets

        if not equal:
            warnings.warn('Mismatch in number of preconditioning tuples.')
            return False

        equal = np.allclose(self.prec, ref.prec, rtol=rtol, atol=atol)

        if not equal:
            warnings.warn('Mismatch in preconditioning parameters.')
            return False

        return True


class Fnetout:
    '''Fortnet xml output instance.'''


    def __init__(self, fname):
        '''Initializes a Precond object.

        Args:

            fname (str): path to working directory to perform regression testing

        '''

        self._fname = fname

        fnetout = ET.parse(self._fname)
        root = fnetout.getroot()
        self.mode = root.find('mode').text.strip()

        if self.mode.startswith(("'", '"')):
            self.mode = self.mode[1:]
        if self.mode.endswith(("'", '"')):
            self.mode = self.mode[:-1]

        if not (self.mode == 'train' or self.mode == 'validate' \
                or self.mode == 'predict'):
            raise FnetoutError('Invalid running mode specification.')

        output = root.find('output')

        self.ndatapoints = int(output.find('ndatapoints').text.strip())
        self.ntargets = int(output.find('ntargets').text.strip())

        if output.find('atomic').text.strip().lower() == 'no':
            self.isatomic = False
        elif output.find('atomic').text.strip().lower() == 'yes':
            self.isatomic = True
        else:
            raise FnetoutError("Invalid specification in 'atomic' tag found.")

        if self.isatomic:
            self.predicts = []
        elif self.mode == 'validate':
            self.predicts = np.empty((self.ndatapoints, 2 * self.ntargets))
        elif self.mode == 'predict':
            self.predicts = np.empty((self.ndatapoints, self.ntargets))

        if self.isatomic:
            for ii in range(self.ndatapoints):
                data = output.find('datapoint' + str(ii + 1))
                natoms = len(list(data))
                if self.mode == 'validate':
                    tmppredicts = np.empty((natoms, 2 * self.ntargets))
                elif self.mode == 'predict':
                    tmppredicts = np.empty((natoms, self.ntargets))
                for jj in range(natoms):
                    tmp = data.find('atom' + str(jj + 1)).text.strip().split()
                    if self.mode == 'validate' and \
                       len(tmp) == 2 * self.ntargets or \
                       (self.mode == 'predict' and \
                        len(tmp) == self.ntargets):
                        tmppredicts[jj, :] = np.array(tmp, dtype=float)
                    else:
                        msg = 'Mismatch in number of given predictions ' + \
                            'and/or targets.'
                        raise FnetoutError(msg)
                self.predicts.append(tmppredicts)
        else:
            for ii in range(self.ndatapoints):
                tmp = output.find(
                    'datapoint' + str(ii + 1)).text.strip().split()
                if self.mode == 'validate' and \
                   len(tmp) == 2 * self.ntargets or \
                   (self.mode == 'predict' and \
                    len(tmp) == self.ntargets):
                    self.predicts[ii, :] = np.array(tmp, dtype=float)
                else:
                    msg = 'Mismatch in number of given predictions ' + \
                        'and/or targets.'
                    raise FnetoutError(msg)


    def equals(self, ref, atol=ATOL, rtol=RTOL):
        '''Checks equality with another reference instance

        Args:

            ref (Acsf): reference instance to compare with
            atol (float): required absolute tolerance
            rtol (float): required relative tolerance

        Returns:

            equal (bool): True, if the two instances are equal within tolerance

        '''

        equal = self.mode == ref.mode

        if not equal:
            warnings.warn('Mismatch in running mode specification.')
            return False

        equal = self.ndatapoints == ref.ndatapoints

        if not equal:
            warnings.warn('Mismatch in total number of datapoints.')
            return False

        equal = self.ntargets == ref.ntargets

        if not equal:
            warnings.warn('Mismatch in number of targets.')
            return False

        equal = self.isatomic == ref.isatomic

        if not equal:
            warnings.warn("Mismatch in 'atomic' tag specification.")
            return False

        if self.isatomic:
            equals = []
            for ii in range(self.ndatapoints):
                equal = np.allclose(self.predicts[ii], ref.predicts[ii],
                                    rtol=rtol, atol=atol)
                equals.append(equal)
            equal = all(equals)
        else:
            equal = np.allclose(self.predicts, ref.predicts, rtol=rtol,
                                atol=atol)

        if not equal:
            warnings.warn('Mismatch in prediction and/or target values.')
            return False

        return True


class TestWithWorkDirError(Exception):
    '''Exception thrown by the TestWithWorkDir class.'''


class AcsfError(Exception):
    '''Exception thrown by the Acsf class.'''


class FnetoutError(Exception):
    '''Exception thrown by the Fnetout class.'''


class NetfileError(Exception):
    '''Exception thrown by the Netfile class.'''


class PrecondError(Exception):
    '''Exception thrown by the Precond class.'''
