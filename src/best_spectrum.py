from typing import Dict

import pandas as pd
from pyteomics import mgf
import spectrum_utils.spectrum as sus


def get_cluster_spectra(mgf_filename: str) -> Dict[int, sus.MsmsSpectrum]:
    """
    Read all spectra from the given MGF file corresponding to a single cluster.

    Parameters
    ----------
    mgf_filename : str
        The file name of the MGF file to be read.

    Returns
    -------
    Dict[int, sus.MsmsSpectrum]
        A dictionary with as keys the scan numbers and as values the
        corresponding spectra.
    """
    return {
        int(spectrum_dict['params']['scans']): sus.MsmsSpectrum(
            spectrum_dict['params']['scans'],
            spectrum_dict['params']['pepmass'][0],
            spectrum_dict['params']['charge'][0],
            spectrum_dict['m/z array'],
            spectrum_dict['intensity array'],
            retention_time=spectrum_dict['params']['rtinseconds'])
        for spectrum_dict in mgf.read(mgf_filename)}


def get_scores(score_filename: str) -> pd.DataFrame:
    """
    Read the PSM scores from the given file.

    Parameters
    ----------
    score_filename : str
        The file name of the MaxQuant identifications file (msms.txt).

    Returns
    -------
    pd.DataFrame
        A `DataFrame` with as columns "Raw file", "Scan number", and "Score".
    """
    return pd.read_csv(score_filename, sep='\t',
                       usecols=['Raw file', 'Scan number', 'Score'])


def get_best_representative(spectra: Dict[int, sus.MsmsSpectrum],
                            scores: pd.DataFrame) -> sus.MsmsSpectrum:
    """
    Get the best representative spectrum with the highest score among the given
    spectra.

    If none of the given spectra have a score a ValueError is raised.

    If two spectra have the same maximum score, the spectrum with the
    alphanumeric first file name with the lowest scan number is selected.

    Parameters
    ----------
    spectra : Dict[int, sus.MsmsSpectrum]
        A dictionary with as keys the scan numbers and as values the
        corresponding spectra in a cluster.
    scores : pd.DataFrame
        A `DataFrame` with the run names in the "Raw file" column, scan number
         in the "Scan number" column, and PSM scores in the "Score" column.

    Returns
    -------
    sus.MsmsSpectrum
        The representative spectrum with the highest score.

    Raises
    ------
    ValueError
        If none of the given spectra have a score.
    """
    scores = (scores[scores['Scan number'].isin(spectra.keys())]
              .sort_values(['Raw file', 'Scan number']))
    if len(scores) == 0:
        raise ValueError('No scores found for the given scan numbers')
    best_scan = scores.iloc[scores['Score'].idxmax()]['Scan number']
    return spectra[best_scan]
