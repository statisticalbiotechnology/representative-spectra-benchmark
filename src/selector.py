import abc
from typing import Dict

import spectrum_utils.spectrum as sus

import ms_io


class RepresentativeSelector(abc.ABC):

    def get_description(self):
        return 'representative'

    def select_representative(
            self, cluster_spectra: Dict[str, sus.MsmsSpectrum]) \
            -> sus.MsmsSpectrum:
        pass


class BestSpectrumRepresentativeSelector(RepresentativeSelector):

    def __init__(self, filename_psms):
        self.psms = ms_io.read_psms(filename_psms)

    def get_description(self):
        return 'best_spectrum'

    def select_representative(
            self, cluster_spectra: Dict[str, sus.MsmsSpectrum]) \
            -> sus.MsmsSpectrum:
        scores = self.psms.loc[self.psms.index.isin(cluster_spectra.keys()),
                               'score']
        if len(scores) == 0:
            raise ValueError('No PSMs found for the given spectra')
        else:
            return cluster_spectra[scores.idxmax()]
