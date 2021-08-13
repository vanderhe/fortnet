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
reference file with the output of the current version of the program.
'''


import os
import warnings
import h5py
import numpy as np


ATOL = 1e-06
RTOL = 1e-05

NETSTAT = 'fortnet.hdf5'
FNETOUT = 'fnetout.hdf5'

ELEMENTSYMBOL = ['h', 'he', 'li', 'be', 'b', 'c', 'n', 'o', 'f', 'ne',
                 'na', 'mg', 'al', 'si', 'p ', 's ', 'cl', 'ar', 'k', 'ca',
                 'sc', 'ti', 'v', 'cr', 'mn', 'fe', 'co', 'ni', 'cu', 'zn',
                 'ga', 'ge', 'as', 'se', 'br', 'kr', 'rb', 'sr', 'y', 'zr',
                 'nb', 'mo', 'tc', 'ru', 'rh', 'pd', 'ag', 'cd', 'in', 'sn',
                 'sb', 'te', 'i', 'xe', 'cs', 'ba', 'la', 'ce', 'pr', 'nd',
                 'pm', 'sm', 'eu', 'gd', 'tb', 'dy', 'ho', 'er', 'tm', 'yb',
                 'lu', 'hf', 'ta', 'w', 're', 'os', 'ir', 'pt', 'au', 'hg',
                 'tl', 'pb', 'bi', 'po', 'at', 'rn', 'fr', 'ra', 'ac', 'th',
                 'pa', 'u', 'np', 'pu', 'am', 'cm', 'bk', 'cf', 'es', 'fm',
                 'md', 'no', 'lr', 'rf', 'db', 'sg', 'bh', 'hs', 'mt', 'ds',
                 'rg', 'cn', 'nh', 'fl', 'mc', 'lv', 'ts', 'og']


class TestWithWorkDir:
    '''Basic Fortnet Regression Test Class.'''


    def __init__(self, wdir):
        '''Initializes a TestWithWorkDir object.

        Args:

            wdir (str): path to working directory to perform regression testing

        '''

        self._wdir = wdir

        # check for the nessecary netstat file
        if not os.path.isfile(os.path.join(self._wdir, '_' + NETSTAT)):
            raise TestWithWorkDirError('Reference netstat file absent.')

        # try to figure out the configuration of the Fortnet run
        self._isfnetout = os.path.isfile(
            os.path.join(self._wdir, '_' + FNETOUT))

        if self._isfnetout and \
           not os.path.isfile(os.path.join(self._wdir, FNETOUT)):
            raise TestWithWorkDirError('Fnetout file absent.')


    def test(self):
        '''Performs regression testing by comparing all relevant files.

        Returns:

            passed (bool): True, if all tests have passed, otherwise False

        '''

        curnetstat = Netfile(os.path.join(self._wdir, NETSTAT))
        refnetstat = Netfile(os.path.join(self._wdir, '_' + NETSTAT))
        netstatpassed = curnetstat.equals(refnetstat)

        if self._isfnetout:
            curfnetout = Fnetout(os.path.join(self._wdir, FNETOUT))
            reffnetout = Fnetout(os.path.join(self._wdir, '_' + FNETOUT))
            fnetoutpassed = curfnetout.equals(reffnetout)
        else:
            fnetoutpassed = True

        passed = netstatpassed and fnetoutpassed

        return passed


class Netfile:
    '''Fortnet network status output instance.'''


    def __init__(self, fname):
        '''Initializes a Netfile object.

        Args:

            fname (str): path to working directory to perform regression testing

        '''

        self._fname = fname
        self.data = dict()
        self.data['bpnn'] = dict()
        self.data['mapping'] = dict()
        self.data['mapping']['preconditioning'] = dict()
        self.data['external'] = dict()

        with h5py.File(self._fname, 'r') as netstatfile:
            netstat = netstatfile['netstat']
            # currently only the BPNN topology is available
            if 'netstat/bpnn' in netstatfile:
                bpnn = netstat['bpnn']
                self.data['nettype'] = 'bpnn'
            else:
                msg = "Error while reading netstat file '" + self._fname + \
                    "'. No network group/information present."
                raise NetfileError(msg)

            self.data['bpnn']['targettype'] = \
                bpnn.attrs.get('targettype').decode('UTF-8').strip()

            if not (self.data['bpnn']['targettype'] == 'atomic' or
                    self.data['bpnn']['targettype'] == 'global'):
                msg = "Error while reading netstat file '" + self._fname + \
                    "'. Unrecognized target type obtained."
                raise NetfileError(msg)
            self.data['bpnn']['atomicnumbers'] = \
                np.sort(np.array(bpnn['atomicnumbers'], dtype=int))

            self.data['bpnn']['nspecies'] = \
                len(self.data['bpnn']['atomicnumbers'])

            for atnum in self.data['bpnn']['atomicnumbers']:
                element = ELEMENTSYMBOL[atnum - 1]
                subnet = bpnn[element + '-subnetwork']
                self.data['bpnn'][element + '-subnetwork'] = dict()
                self.data['bpnn'][element + '-subnetwork']['topology'] = \
                    np.array(subnet['topology'], dtype=int)
                self.data['bpnn'][element + '-subnetwork']['nlayer'] = \
                    len(self.data['bpnn'][element + '-subnetwork']['topology'])
                self.data['bpnn'][element + '-subnetwork']['activation'] = \
                    subnet.attrs.get('activation').decode('UTF-8').strip()
                self.data['bpnn'][element + '-subnetwork']['element'] = \
                    subnet.attrs.get('element')[0]
                self.data['bpnn'][element + '-subnetwork']['nlayer'] = \
                    subnet.attrs.get('nlayer')[0]
                for ilayer in range(self.data['bpnn'] \
                                    [element + '-subnetwork']['nlayer'] - 1):
                    layername = 'layer' + str(ilayer + 1)
                    layer = subnet[layername]
                    self.data['bpnn'][element + '-subnetwork'][layername] \
                        = dict()
                    self.data['bpnn'][element + '-subnetwork'][layername] \
                        ['weights'] = np.array(layer['weights'], dtype=float)
                    self.data['bpnn'][element + '-subnetwork'][layername] \
                        ['bias'] = np.array(layer['bias'], dtype=float)

            # inquire structural ACSF mappings
            if 'netstat/mapping' in netstatfile:
                mapping = netstat['mapping']
                self.data['mapping']['type'] = \
                    mapping.attrs.get('type').decode('UTF-8').strip()
                if not self.data['mapping']['type'] == 'acsf':
                    msg = "Error while reading netstat file '" + self._fname + \
                        "'. Unrecognized mapping type obtained."
                    raise NetfileError(msg)
            else:
                self.data['mapping']['type'] = None
                self.data['mapping']['nfunctions'] = 0
                self.data['mapping']['preconditioning']['type'] = None

            if self.data['mapping']['type'] is not None:
                self.data['mapping']['nfunctions'] = \
                    mapping.attrs.get('nfunctions')
                if len(self.data['mapping']['nfunctions']) == 1:
                    self.data['mapping']['nfunctions'] = \
                        self.data['mapping']['nfunctions'][0]
                else:
                    msg = "Error while reading netstat file '" + self._fname + \
                        "'. Unrecognized number of ACSF functions obtained."
                    raise NetfileError(msg)
                if 'preconditioning' in mapping:
                    precond = mapping['preconditioning']
                    self.data['mapping']['preconditioning']['type'] = \
                        precond.attrs.get('type').decode('UTF-8').strip()
                    if not self.data['mapping']['preconditioning']['type'] \
                       == 'zscore':
                        msg = "Error while reading netstat file '" + \
                            self._fname + \
                            "'. Unrecognized mapping type obtained."
                        raise NetfileError(msg)
                else:
                    self.data['mapping']['preconditioning']['type'] = None

            if self.data['mapping']['type'] is not None:
                for ifunc in range(self.data['mapping']['nfunctions']):
                    funcname = 'function' + str(ifunc + 1)
                    function = mapping[funcname]
                    self.data['mapping'][funcname] = dict()
                    self.data['mapping'][funcname]['type'] = \
                        function.attrs.get('type').decode('UTF-8').strip()
                    self.data['mapping'][funcname]['atomicnumbers'] = \
                        np.array(function.attrs.get('atomicnumbers'), dtype=int)
                    self.data['mapping'][funcname] \
                        ['atomid'] = function.attrs.get('atomid')[0]
                    self.data['mapping'][funcname] \
                        ['cutoff'] = function.attrs.get('cutoff')[0]
                    if 'eta' in function.attrs:
                        self.data['mapping'][funcname] \
                            ['eta'] = function.attrs.get('eta')[0]
                    else:
                        self.data['mapping'][funcname] \
                            ['eta'] = None
                    if 'lambda' in function.attrs:
                        self.data['mapping'][funcname] \
                            ['lambda'] = function.attrs.get('lambda')[0]
                    else:
                        self.data['mapping'][funcname] \
                            ['lambda'] = None
                    if 'xi' in function.attrs:
                        self.data['mapping'][funcname] \
                            ['xi'] = function.attrs.get('xi')[0]
                    else:
                        self.data['mapping'][funcname] \
                            ['xi'] = None
                    if 'kappa' in function.attrs:
                        self.data['mapping'][funcname] \
                            ['kappa'] = function.attrs.get('kappa')[0]
                    else:
                        self.data['mapping'][funcname] \
                            ['kappa'] = None
                    tmp = function.attrs.get('tradial')[0]
                    if tmp == 1:
                        self.data['mapping'][funcname]['tradial'] = True
                    elif tmp == 0:
                        self.data['mapping'][funcname]['tradial'] = False
                    else:
                        msg = "Error while reading netstat file '" + \
                            self._fname + \
                            "'. Invalid mapping value tRadial obtained."
                        raise NetfileError(msg)
                    tmp = function.attrs.get('tangular')[0]
                    if tmp == 1:
                        self.data['mapping'][funcname]['tangular'] = True
                    elif tmp == 0:
                        self.data['mapping'][funcname]['tangular'] = False
                    else:
                        msg = "Error while reading netstat file '" + \
                            self._fname + \
                            "'. Invalid mapping value tAngular obtained."
                        raise NetfileError(msg)

            if self.data['mapping']['preconditioning']['type'] is not None:
                self.data['mapping']['preconditioning']['means'] = \
                    np.array(precond['means'], dtype=float)
                self.data['mapping']['preconditioning']['variances'] = \
                    np.array(precond['variances'], dtype=float)

            # inquire external atomic input features
            if 'netstat/external' in netstatfile:
                external = netstat['external']
                self.data['external']['nextfeatures'] = \
                    external.attrs.get('nextfeatures')
                if len(self.data['external']['nextfeatures']) == 1:
                    self.data['external']['nextfeatures'] = \
                        self.data['external']['nextfeatures'][0]
                else:
                    msg = "Error while reading netstat file '" + self._fname + \
                        "'. Unrecognized number of external features obtained."
                    raise NetfileError(msg)
                self.data['external']['indices'] = \
                    np.array(external['indices'], dtype=int)
            else:
                self.data['external']['nextfeatures'] = 0
                self.data['external']['indices'] = None


    def equals(self, ref, atol=ATOL, rtol=RTOL):
        '''Checks equality with another reference instance

        Args:

            ref (Acsf): reference instance to compare with
            atol (float): required absolute tolerance
            rtol (float): required relative tolerance

        Returns:

            equal (bool): True, if the two instances are equal within tolerance

        '''

        equal = self.data['nettype'] == ref.data['nettype']

        if not equal:
            warnings.warn('Mismatch in network type.')
            return False

        equal = self.data['bpnn']['targettype'] == \
            ref.data['bpnn']['targettype']

        if not equal:
            warnings.warn('Mismatch in target type.')
            return False

        equal = self.data['bpnn']['nspecies'] == ref.data['bpnn']['nspecies']

        if not equal:
            warnings.warn('Mismatch in number of BPNN species.')
            return False

        equal = np.allclose(self.data['bpnn']['atomicnumbers'],
                            ref.data['bpnn']['atomicnumbers'],
                            rtol=rtol, atol=atol)

        if not equal:
            warnings.warn('Mismatch in atomic numbers of sub-networks.')
            return False

        for atnum in self.data['bpnn']['atomicnumbers']:
            element = ELEMENTSYMBOL[atnum - 1]
            equal = \
                (self.data['bpnn'][element + '-subnetwork']['topology'] \
                 == ref.data['bpnn'][element + '-subnetwork']['topology']) \
                 .all()
            if not equal:
                warnings.warn('Mismatch in topologies of ' + element + \
                              '-subnetwork.')
                return False
            equal = (self.data['bpnn'][element + '-subnetwork']['nlayer'] == \
                     ref.data['bpnn'][element + '-subnetwork']['nlayer'])
            if not equal:
                warnings.warn('Mismatch in number of layers in the ' + element \
                              + '-subnetwork.')
                return False
            equal = (self.data['bpnn'][element + '-subnetwork']['element'] == \
                     ref.data['bpnn'][element + '-subnetwork']['element'])
            if not equal:
                warnings.warn('Mismatch in element specification in the ' + \
                              element + '-subnetwork.')
                return False
            equal = \
                (self.data['bpnn'][element + '-subnetwork']['activation'] == \
                 ref.data['bpnn'][element + '-subnetwork']['activation'])
            if not equal:
                warnings.warn('Mismatch in activation function in the ' + \
                              element + '-subnetwork.')
                return False
            for ilayer in range(self.data['bpnn'] \
                                [element + '-subnetwork']['nlayer'] - 1):
                equal = np.allclose(
                    self.data['bpnn'][element + '-subnetwork'] \
                    ['layer' + str(ilayer + 1)]['weights'],
                    ref.data['bpnn'][element + '-subnetwork'] \
                    ['layer' + str(ilayer + 1)]['weights'],
                    rtol=rtol, atol=atol)
                if not equal:
                    warnings.warn('Mismatch in weight parameters of ' + \
                                  element + '-subnetwork (layer ' \
                                  + str(ilayer + 1) + ').')
                    return False
                equal = np.allclose(
                    self.data['bpnn'][element + '-subnetwork'] \
                    ['layer' + str(ilayer + 1)]['bias'],
                    ref.data['bpnn'][element + '-subnetwork']['layer' + \
                    str(ilayer + 1)]['bias'],
                    rtol=rtol, atol=atol)
                if not equal:
                    warnings.warn('Mismatch in bias parameters of ' + \
                                  element + '-subnetwork (layer ' \
                                  + str(ilayer + 1) + ').')
                    return False

        if self.data['mapping']['type'] is not None:
            equal = (self.data['mapping']['type'] == \
                     ref.data['mapping']['type'])

            equal = (self.data['mapping']['nfunctions'] == \
                     ref.data['mapping']['nfunctions'])

            for ifunc in range(self.data['mapping']['nfunctions']):
                funcname = 'function' + str(ifunc + 1)

                equal = (self.data['mapping'][funcname]['type'] == \
                         ref.data['mapping'][funcname]['type'])
                if not equal:
                    warnings.warn('Mismatch in G-function type of function ' + \
                                  str(ifunc + 1) + '.')
                    return False

                equal = np.allclose(
                    self.data['mapping'][funcname]['atomicnumbers'],
                    ref.data['mapping'][funcname]['atomicnumbers'],
                    rtol=rtol, atol=atol)
                if not equal:
                    warnings.warn('Mismatch in atomic numbers of function ' + \
                                  str(ifunc + 1) + '.')
                    return False

                equal = np.allclose(
                    [self.data['mapping'][funcname]['atomid']],
                    [ref.data['mapping'][funcname]['atomid']],
                    rtol=rtol, atol=atol)
                if not equal:
                    warnings.warn('Mismatch in atomId of function ' + \
                                  str(ifunc + 1) + '.')
                    return False

                equal = np.allclose(
                    [self.data['mapping'][funcname]['cutoff']],
                    [ref.data['mapping'][funcname]['cutoff']],
                    rtol=rtol, atol=atol)
                if not equal:
                    warnings.warn('Mismatch in cutoff-parameter of function ' \
                                  + str(ifunc + 1) + '.')
                    return False

                if self.data['mapping'][funcname]['eta'] is not None:
                    equal = np.allclose(
                        [self.data['mapping'][funcname]['eta']],
                        [ref.data['mapping'][funcname]['eta']],
                        rtol=rtol, atol=atol)
                else:
                    equal = True
                if not equal:
                    warnings.warn('Mismatch in eta-parameter of function ' + \
                                  str(ifunc + 1) + '.')
                    return False

                if self.data['mapping'][funcname]['lambda'] is not None:
                    equal = np.allclose(
                        [self.data['mapping'][funcname]['lambda']],
                        [ref.data['mapping'][funcname]['lambda']],
                        rtol=rtol, atol=atol)
                else:
                    equal = True
                if not equal:
                    warnings.warn('Mismatch in lambda-parameter of function ' \
                                  + str(ifunc + 1) + '.')
                    return False

                if self.data['mapping'][funcname]['xi'] is not None:
                    equal = np.allclose(
                        [self.data['mapping'][funcname]['xi']],
                        [ref.data['mapping'][funcname]['xi']],
                        rtol=rtol, atol=atol)
                else:
                    equal = True
                if not equal:
                    warnings.warn('Mismatch in xi-parameter of function ' + \
                                  str(ifunc + 1) + '.')
                    return False

                if self.data['mapping'][funcname]['kappa'] is not None:
                    equal = np.allclose(
                        [self.data['mapping'][funcname]['kappa']],
                        [ref.data['mapping'][funcname]['kappa']],
                        rtol=rtol, atol=atol)
                else:
                    equal = True
                if not equal:
                    warnings.warn('Mismatch in kappa-parameter of function ' + \
                                  str(ifunc + 1) + '.')
                    return False

                equal = (self.data['mapping'][funcname]['tradial'] == \
                         ref.data['mapping'][funcname]['tradial'])
                if not equal:
                    warnings.warn('Mismatch in tRadial specification of ' + \
                                  'function ' + str(ifunc + 1) + '.')
                    return False

                equal = (self.data['mapping'][funcname]['tangular'] == \
                         ref.data['mapping'][funcname]['tangular'])
                if not equal:
                    warnings.warn('Mismatch in tAngular specification of ' + \
                                  'function ' + str(ifunc + 1) + '.')
                    return False

        if self.data['mapping']['preconditioning']['type'] is not None:
            equal = (self.data['mapping']['preconditioning']['type'] == \
                     ref.data['mapping']['preconditioning']['type'])
            if not equal:
                warnings.warn('Mismatch in ACSF preconditioning type.')
                return False

            equal = np.allclose(
                self.data['mapping']['preconditioning']['means'],
                ref.data['mapping']['preconditioning']['means'],
                rtol=rtol, atol=atol)
            if not equal:
                warnings.warn('Mismatch z-score mean values.')
                return False

            equal = np.allclose(
                self.data['mapping']['preconditioning']['variances'],
                ref.data['mapping']['preconditioning']['variances'],
                rtol=rtol, atol=atol)
            if not equal:
                warnings.warn('Mismatch in z-score variances.')
                return False

        if self.data['external']['indices'] is not None:
            equal = (self.data['external']['nextfeatures'] == \
                     ref.data['external']['nextfeatures'])
            if not equal:
                warnings.warn('Mismatch in number of external features.')
                return False

            equal = (self.data['external']['indices'] == \
                     ref.data['external']['indices']).all()
            if not equal:
                warnings.warn('Mismatch in external feature indices.')
                return False

        return True


class Fnetout:
    '''Fortnet output instance.'''


    def __init__(self, fname):
        '''Initializes a Fnetout object.

        Args:

            fname (str): path to working directory to perform regression testing

        '''

        self._fname = fname
        self.data = dict()
        self.data['fnetout'] = dict()
        self.data['fnetout']['output'] = dict()

        with h5py.File(self._fname, 'r') as fnetoutfile:
            fnetout = fnetoutfile['fnetout']
            self.data['fnetout']['mode'] = \
                fnetout.attrs.get('mode').decode('UTF-8').strip()
            if not (self.data['fnetout']['mode'] == 'validate' or \
                    self.data['fnetout']['mode'] == 'predict'):
                raise FnetoutError('Invalid running mode specification.')
            output = fnetout['output']
            self.data['fnetout']['output']['ndatapoints'] = \
                output.attrs.get('ndatapoints')
            if len(self.data['fnetout']['output']['ndatapoints']) == 1:
                self.data['fnetout']['output']['ndatapoints'] = \
                    self.data['fnetout']['output']['ndatapoints'][0]
            else:
                msg = "Error while reading fnetout file '" + self._fname + \
                    "'. Unrecognized number of datapoints obtained."
                raise FnetoutError(msg)
            self.data['fnetout']['output']['targettype'] = \
                output.attrs.get('targettype').decode('UTF-8').strip()
            if not (self.data['fnetout']['output']['targettype'] == 'atomic' \
                    or self.data['fnetout']['output']['targettype'] \
                    == 'global'):
                raise FnetoutError('Invalid running mode specification.')

            for idata in range(self.data['fnetout']['output']['ndatapoints']):
                dataname = 'datapoint' + str(idata + 1)
                self.data['fnetout']['output'][dataname] = dict()
                if self.data['fnetout']['mode'] == 'validate':
                    ntargets = int(
                        np.shape(np.array(output[dataname]['output']))[1] / 2)
                    self.data['fnetout']['output'][dataname]['predictions'] = \
                        np.array(output[dataname]['output'],
                                 dtype=float)[:ntargets, :]
                    self.data['fnetout']['output'][dataname]['targets'] = \
                        np.array(output[dataname]['output'],
                                 dtype=float)[ntargets:, :]
                else:
                    self.data['fnetout']['output'][dataname]['predictions'] = \
                        np.array(output[dataname]['output'], dtype=float)
                    self.data['fnetout']['output'][dataname]['targets'] = None


    def equals(self, ref, atol=ATOL, rtol=RTOL):
        '''Checks equality with another reference instance

        Args:

            ref (Acsf): reference instance to compare with
            atol (float): required absolute tolerance
            rtol (float): required relative tolerance

        Returns:

            equal (bool): True, if the two instances are equal within tolerance

        '''

        equal = self.data['fnetout']['mode'] == ref.data['fnetout']['mode']

        if not equal:
            warnings.warn('Mismatch in running mode specification.')
            return False

        equal = self.data['fnetout']['output']['ndatapoints'] == \
            ref.data['fnetout']['output']['ndatapoints']

        if not equal:
            warnings.warn('Mismatch in total number of datapoints.')
            return False

        equal = self.data['fnetout']['output']['targettype'] == \
            ref.data['fnetout']['output']['targettype']

        if not equal:
            warnings.warn('Mismatch in target type specification.')
            return False

        for idata in range(self.data['fnetout']['output']['ndatapoints']):
            dataname = 'datapoint' + str(idata + 1)
            equal = np.allclose(
                self.data['fnetout']['output'][dataname]['predictions'],
                ref.data['fnetout']['output'][dataname]['predictions'],
                rtol=rtol, atol=atol)
            if not equal:
                warnings.warn('Mismatch in predictions of datapoint ' \
                              + str(idata + 1) + '.')
                return False
            if self.data['fnetout']['mode'] == 'validate':
                equal = np.allclose(
                    self.data['fnetout']['output'][dataname]['targets'],
                    ref.data['fnetout']['output'][dataname]['targets'],
                    rtol=rtol, atol=atol)
                if not equal:
                    warnings.warn('Mismatch in targets of datapoint ' \
                                  + str(idata + 1) + '.')
                    return False

        return True


class TestWithWorkDirError(Exception):
    '''Exception thrown by the TestWithWorkDir class.'''


class NetfileError(Exception):
    '''Exception thrown by the Netfile class.'''


class FnetoutError(Exception):
    '''Exception thrown by the Fnetout class.'''
