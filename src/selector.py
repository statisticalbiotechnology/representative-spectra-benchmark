import abc
import copy
import functools
import math
from typing import Dict

import numpy as np
import spectrum_utils.spectrum as sus

import metrics
import ms_io


class RepresentativeSelector(abc.ABC):
    """
    Abstract parent class for strategies to select representative spectra.
    """

    def get_description(self) -> str:
        """
        Brief description of the selection strategy.

        Returns
        -------
        str
            General description: 'representative'.
        """
        return 'representative'

    def select_representative(
            self, cluster_spectra: Dict[str, sus.MsmsSpectrum]) \
            -> sus.MsmsSpectrum:
        """
        Select a representative spectrum given a collection of spectra.

        Parameters
        ----------
        cluster_spectra : Dict[str, sus.MsmsSpectrum]
            A collection of spectra forming a cluster. The dictionary consists
            of spectrum key -> MsmsSpectrum object pairs.

        Returns
        -------
        sus.MsmsSpectrum
            The representative spectrum for the given collection of spectra.
        """
        pass


class BestSpectrumRepresentativeSelector(RepresentativeSelector):
    """
    Select as cluster representative the spectrum with the highest search
    engine score.
    """

    def __init__(self, filename_psms: str, higher_is_better: bool = True):
        """
        Initialize the selector with the search engine scores for all spectra
        to be considered.

        Parameters
        ----------
        filename_psms : str
            The file name of the search engine PSMs. Identifications for all
            spectra to be successfully considered need to be included in this
            file.
        """
        self.psms = ms_io.read_psms(filename_psms)
        self.higher_is_better = higher_is_better

    def get_description(self):
        """
        Brief description of the selection strategy.

        Returns
        -------
        str
            Description: 'best_spectrum'.
        """
        return 'best_spectrum'

    def select_representative(
            self, cluster_spectra: Dict[str, sus.MsmsSpectrum]) \
            -> sus.MsmsSpectrum:
        """
        Select the representative spectrum with the highest search engine score
        in the given a collection of spectra.

        Parameters
        ----------
        cluster_spectra : Dict[str, sus.MsmsSpectrum]
            A collection of spectra forming a cluster. The dictionary consists
            of spectrum key -> MsmsSpectrum object pairs.

        Returns
        -------
        sus.MsmsSpectrum
            The representative spectrum with the highest search engine score
            for the given collection of spectra.
        """
        scores = self.psms.loc[self.psms.index.isin(cluster_spectra.keys()),
                               'score']
        if len(scores) == 0:
            raise ValueError('No PSMs found for the given spectra')
        else:
            return cluster_spectra[scores.idxmax() if self.higher_is_better
                                   else scores.idxmin()]


class MostSimilarRepresentativeSelector(RepresentativeSelector):
    """
    Select as cluster representative the most similar spectrum that has the
    maximum average similarity to all other spectra in the cluster.
    """

    def __init__(self, sim: str, fragment_mz_tolerance: float):
        """
        Initialize the selector with the given similarity measure.

        Supported similarity measures:
        - 'dot'

        Parameters
        ----------
        sim : str
            The similarity measure used to compare spectra with each other
            (default: dot product).
        fragment_mz_tolerance : float
            The fragment m/z tolerance used during spectrum comparisons.
        """
        self.sim = sim
        if self.sim == 'dot':
            self.compare_spectra = functools.partial(
                metrics.dot, fragment_mz_tolerance=fragment_mz_tolerance)
        else:
            raise ValueError('Unknown spectrum similarity method')

    def get_description(self):
        """
        Brief description of the selection strategy.

        Returns
        -------
        str
            Description: 'most_similar_' + similarity measure.
        """
        return f'most_similar_{self.sim}'

    def select_representative(
            self, cluster_spectra: Dict[str, sus.MsmsSpectrum]) \
            -> sus.MsmsSpectrum:
        """
        Select the representative spectrum with the maximum average similarity
        to all other spectra in the cluster.

        Parameters
        ----------
        cluster_spectra : Dict[str, sus.MsmsSpectrum]
            A collection of spectra forming a cluster. The dictionary consists
            of spectrum key -> MsmsSpectrum object pairs.

        Returns
        -------
        sus.MsmsSpectrum
            The representative spectrum with the maximum average similarity to
            the other spectra in the given collection of spectra.
        """
        spectra_keys = list(cluster_spectra.keys())
        if len(spectra_keys) == 1:
            return cluster_spectra[spectra_keys[0]]
        # Compute pairwise similarities.
        sim_matrix = np.zeros((len(spectra_keys), len(spectra_keys)))
        for i in range(len(spectra_keys)):
            for j in range(i, len(spectra_keys)):
                sim_matrix[i, j] = sim_matrix[j, i] = self.compare_spectra(
                    cluster_spectra[spectra_keys[i]],
                    cluster_spectra[spectra_keys[j]])
        # Find the spectrum with the maximum similarity to all other spectra.
        return cluster_spectra[spectra_keys[sim_matrix.sum(axis=0).argmax()]]


class BinningRepresentativeSelector(RepresentativeSelector):
    """
    Select as cluster representative a binned average of all the spectra in the
    cluster.
    """

    def __init__(self, min_mz: float, max_mz: float, bin_size: float,
                 peak_quorum: float, edge_case_threshold: float):
        """
        Initialize the selector with the information needed to compute average
        binned representatives.

        Parameters
        ----------
        min_mz : float
            The minimum m/z to consider.
        max_mz : float
            The maximum m/z to consider.
        bin_size : float
            Bin size in m/z.
        peak_quorum : float
            Relative number of spectra in a cluster that need to contain a peak
            for it to be included in the representative spectrum.
        edge_case_threshold : float
            Try to correct m/z edge cases where the m/z is closer to the bin
            edge than the given relative bin size threshold.
        """
        self.min_mz = min_mz
        self.max_mz = max_mz
        self.bin_size = bin_size
        self.peak_quorum = peak_quorum
        self.edge_case_threshold = edge_case_threshold

    def get_description(self):
        """
        Brief description of the selection strategy.

        Returns
        -------
        str
            Description: 'bin_min{min_mz}_max{max_mz}_bin{bin_size}_
                          quorum{peak_quorum}_edge{edge_case_threshold}'.
        """
        return (f'bin_min{self.min_mz}_max{self.max_mz}_bin{self.bin_size}_'
                f'quorum{self.peak_quorum}_edge{self.edge_case_threshold}'
                .replace('.', ''))

    def select_representative(
            self, cluster_spectra: Dict[str, sus.MsmsSpectrum]) \
            -> sus.MsmsSpectrum:
        """
        Compile a representative spectrum that is the binned average of the
        given spectra.

        Parameters
        ----------
        cluster_spectra : Dict[str, sus.MsmsSpectrum]
            A collection of spectra forming a cluster. The dictionary consists
            of spectrum key -> MsmsSpectrum object pairs.

        Returns
        -------
        sus.MsmsSpectrum
            A binned average representative spectrum from the collection of
            spectra.
        """
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
