import sys
from typing import Iterable

import numba as nb
import numpy as np
import spectrum_utils.spectrum as sus
from pyteomics import parser


mz_unit = 1.000508  # Space in Th between fragments
mz_space = mz_unit * .005  # Resolution approximately 0.005 Da


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


def fraction_of_by(representative_spectrum, cluster_members=[]):
    if not representative_spectrum.peptide:
        return 0.0
    fragment_tol_mass = 0.005
    fragment_tol_mode = 'Da'
    spectrum = (representative_spectrum.remove_precursor_peak(fragment_tol_mass, fragment_tol_mode)
                .annotate_peptide_fragments(fragment_tol_mass, fragment_tol_mode,
                                            ion_types='by'))
    current, by_current = 0., 0.
    for ix in range(len(spectrum.intensity)):
        current += spectrum.intensity[ix]
        if spectrum.annotation[ix] is not None:
            by_current += spectrum.intensity[ix]
    if current > 0.:
        return by_current / current
    else:
        return 0.0


def fraction_of_by_seq(peptide_seq, precursor_mz, precursor_charge, mz, intensity):
    if not parser.fast_valid(peptide_seq):
        print("Invalid peptide sequence encountered", file=sys.stderr)
        return 0.0
    spec = sus.MsmsSpectrum(
        peptide_seq, precursor_mz=precursor_mz, precursor_charge=precursor_charge, mz=mz, intensity=intensity,
        peptide=peptide_seq)
    return fraction_of_by(spec)
