import abc
import copy
import math
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


class BinningRepresentativeSelector(RepresentativeSelector):

    def __init__(self, min_mz: float = 100., max_mz: float = 2000.,
                 bin_size: float = 0.02, peak_quorum: float = 0.25,
                 edge_case_threshold: float = 0.5):
        self.min_mz = min_mz
        self.max_mz = max_mz
        self.bin_size = bin_size
        self.peak_quorum = peak_quorum
        self.edge_case_threshold = edge_case_threshold

    def get_description(self):
        return 'bin'

    def select_representative(
            self, cluster_spectra: Dict[str, sus.MsmsSpectrum]) \
            -> sus.MsmsSpectrum:

        num_bins = math.ceil((self.max_mz - self.min_mz) / self.bin_size)
        mzs = np.zeros(num_bins, dtype=np.float32)
        intensities = np.zeros(num_bins, dtype=np.float32)
        n_peaks = np.zeros(num_bins, dtype=np.uint32)
        precursor_mzs, precursor_charges = [], []
        # Aggregate spectra peaks.
        for spectrum in cluster_spectra.values():
            spectrum = copy.copy(spectrum).set_mz_range(
                self.min_mz, self.max_mz)
            bin_array = np.floor((spectrum.mz - self.min_mz)
                                 / self.bin_size).astype(np.uint32)
            n_peaks[bin_array] += 1
            intensities[bin_array] += spectrum.intensity
            mzs[bin_array] += spectrum.mz
            precursor_mzs.append(spectrum.precursor_mz)
            precursor_charges.append(spectrum.precursor_charge)

        # Verify that all precursor charges are the same.
        if not all(charge == precursor_charges[0]
                   for charge in precursor_charges):
            raise ValueError('Spectra in a cluster have different precursor '
                             'charges')

        # Try to handle the case where a single peak is split on a bin
        # boundary.
        mz_temp = np.copy(mzs)
        mz_temp_mask = n_peaks > 0
        mz_temp[mz_temp_mask] /= n_peaks[mz_temp_mask]
        # Subtract the mzs from their previous mz.
        mz_delta = np.diff(mz_temp)
        mz_delta[-1] = 0
        # Handle cases where the deltas are smaller than the thresholded bin
        # size.
        mz_delta_small_index = np.nonzero(
            (mz_delta > 0) &
            (mz_delta < self.bin_size * self.edge_case_threshold))[0]
        if len(mz_delta_small_index) > 0:
            # Consolidate all the split mzs, intensities, n_peaks into one bin.
            mzs[mz_delta_small_index] += mzs[mz_delta_small_index + 1]
            mzs[mz_delta_small_index + 1] = 0
            intensities[mz_delta_small_index] += \
                intensities[mz_delta_small_index + 1]
            intensities[mz_delta_small_index + 1] = 0
            n_peaks[mz_delta_small_index] += n_peaks[mz_delta_small_index + 1]
            n_peaks[mz_delta_small_index + 1] = 0

        # Determine how many peaks need to be present to keep a final peak.
        peak_quorum_int = math.ceil(len(cluster_spectra) * self.peak_quorum)
        mask = n_peaks >= peak_quorum_int
        # Take the mean of all peaks per bin.
        mzs = mzs[mask] / n_peaks[mask]
        intensities = intensities[mask] / n_peaks[mask]
        precursor_mz = np.mean(precursor_mzs)
        precursor_charge = precursor_charges[0]

        return sus.MsmsSpectrum('binned_representative', precursor_mz,
                                precursor_charge, mzs, intensities)
