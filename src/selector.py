import abc
import copy
from typing import Dict

import numpy as np
import pyopenms
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


def dist(spectrum1: sus.MsmsSpectrum, spectrum2: sus.MsmsSpectrum,
         method: str = 'dot'):
    if method == 'dot':
        spectrum1 = copy.copy(spectrum1).scale_intensity(max_intensity=1)
        spec1 = pyopenms.MSSpectrum()
        spec1.set_peaks([spectrum1.mz, spectrum1.intensity])
        spectrum2 = copy.copy(spectrum2).scale_intensity(max_intensity=1)
        spec2 = pyopenms.MSSpectrum()
        spec2.set_peaks([spectrum2.mz, spectrum2.intensity])
        # Third parameter of xCorrelationPrescore is bin size in Da.
        xcorr = pyopenms.XQuestScores().xCorrelationPrescore(spec1, spec2,
                                                             0.005)
        return 1 - xcorr
    else:
        return 0


class MostSimilarRepresentativeSelector(RepresentativeSelector):

    def __init__(self, distance: str = 'dot'):
        self.distance = distance if distance is not None else 'dot'

    def get_description(self):
        return f'most_similar_{self.distance}'

    def select_representative(
            self, cluster_spectra: Dict[str, sus.MsmsSpectrum]) \
            -> sus.MsmsSpectrum:
        spectra_keys = list(cluster_spectra.keys())
        if len(spectra_keys) == 1:
            return cluster_spectra[spectra_keys[0]]
        # Compute pairwise distances.
        dist_matrix = np.zeros((len(spectra_keys), len(spectra_keys)))
        for i in range(len(spectra_keys)):
            for j in range(i, len(spectra_keys)):
                dist_matrix[i, j] = dist_matrix[j, i] = dist(
                    cluster_spectra[spectra_keys[i]],
                    cluster_spectra[spectra_keys[j]],
                    self.distance)
        # Find the spectrum with the minimal distance to all other spectra.
        return cluster_spectra[spectra_keys[dist_matrix.sum(axis=0).argmin()]]
