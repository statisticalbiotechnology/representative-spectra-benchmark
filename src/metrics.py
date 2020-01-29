import copy
from typing import Iterable

import numba as nb
import numpy as np
import spectrum_utils.spectrum as sus


def dot(spectrum1: sus.MsmsSpectrum, spectrum2: sus.MsmsSpectrum,
        fragment_mz_tolerance: float) -> float:
    """
    Compute the dot product between the given spectra.

    Parameters
    ----------
    spectrum1 : sus.MsmsSpectrum
        The first spectrum.
    spectrum2 : sus.MsmsSpectrum
        The second spectrum.
    fragment_mz_tolerance : float
        The fragment m/z tolerance used to match peaks.

    Returns
    -------
    float
        The dot product similarity between the given spectra.
    """
    return _dot(spectrum1.mz, _norm_intensity(np.copy(spectrum1.intensity)),
                spectrum2.mz, _norm_intensity(np.copy(spectrum2.intensity)),
                fragment_mz_tolerance)


@nb.njit
def _norm_intensity(spectrum_intensity: np.ndarray) -> np.ndarray:
    """
    Normalize spectrum peak intensities.

    Parameters
    ----------
    spectrum_intensity : np.ndarray
        The spectrum peak intensities to be normalized.

    Returns
    -------
    np.ndarray
        The normalized peak intensities.
    """
    return spectrum_intensity / np.linalg.norm(spectrum_intensity)


@nb.njit
def _dot(mz: np.ndarray, intensity: np.ndarray, mz_other: np.ndarray,
         intensity_other: np.ndarray, fragment_mz_tol: float) -> float:
    """
    Compute the dot product between the given spectra.

    Note: Spectrum intensities should be normalized prior to computing the dot
    product.

    Parameters
    ----------
    mz : np.ndarray
        The first spectrum's m/z values.
    intensity : np.ndarray
        The first spectrum's intensity values.
    mz_other : np.ndarray
        The second spectrum's m/z values.
    intensity_other : np.ndarray
        The second spectrum's intensity values.
    fragment_mz_tol : float
        The fragment m/z tolerance used to match peaks in both spectra with
        each other.

    Returns
    -------
    float
        The dot product between both spectra.
    """
    fragment_i, fragment_other_i, score = 0, 0, 0.
    for fragment_i in range(len(mz)):
        while (fragment_other_i < len(mz_other) - 1 and
               mz_other[fragment_other_i] < mz[fragment_i] - fragment_mz_tol):
            fragment_other_i += 1
        if (abs(mz[fragment_i] - mz_other[fragment_other_i]) <= fragment_mz_tol
                and fragment_other_i < len(mz_other)):
            score += intensity[fragment_i] * intensity_other[fragment_other_i]
            fragment_other_i += 1
    return score


def avg_dot(representative: sus.MsmsSpectrum,
            cluster_spectra: Iterable[sus.MsmsSpectrum],
            fragment_mz_tolerance: float) -> float:
    """
    Compute the average dot product between the cluster representative and all
    cluster members.

    Parameters
    ----------
    representative : sus.MsmsSpectrum
        The cluster representative spectrum.
    cluster_spectra : Iterable[sus.MsmsSpectrum]
        The cluster member spectra.
    fragment_mz_tolerance : float
        Fragment m/z tolerance used during spectrum comparison.

    Returns
    -------
    float
        The average dot product between the cluster representative and all
        cluster members.
    """
    return np.mean([dot(representative, spectrum, fragment_mz_tolerance)
                    for spectrum in cluster_spectra])


def fraction_by(representative: sus.MsmsSpectrum,
                cluster_spectra: Iterable[sus.MsmsSpectrum],
                fragment_mz_tolerance: float) -> float:
    """
    Compute the fraction of intensity that is explained by the b and y-ions of
    the representative spectrum.

    This will be 0 if no peptide sequence is associated with the representative
    spectrum.

    Parameters
    ----------
    representative : sus.MsmsSpectrum
        The cluster representative spectrum.
    cluster_spectra : Iterable[sus.MsmsSpectrum]
        The cluster member spectra. Ignored.
    fragment_mz_tolerance : float
        Fragment m/z tolerance used to annotate the peaks of the representative
        spectrum.

    Returns
    -------
    float
        The fraction of intensity that is explained by the b and y-ions of the
        representative spectrum.
    """
    if representative.peptide is None:
        return 0
    representative = (copy.copy(representative)
                      .remove_precursor_peak(fragment_mz_tolerance, 'Da')
                      .annotate_peptide_fragments(fragment_mz_tolerance, 'Da'))
    annotated_peaks = [i for i, annot in enumerate(representative.annotation)
                       if annot is not None]
    return (representative.intensity[annotated_peaks].sum()
            / representative.intensity.sum())
