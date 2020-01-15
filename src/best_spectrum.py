import collections
import sys
from typing import Dict, Iterable, List

import pandas as pd
from pyteomics import mgf
import spectrum_utils.spectrum as sus


def get_cluster_spectra(mgf_filename: str) -> Dict[str, sus.MsmsSpectrum]:
    """
    Read all spectra from the given MGF file corresponding to a single cluster.

    Parameters
    ----------
    mgf_filename : str
        The file name of the MGF file to be read.

    Returns
    -------
    Dict[str, sus.MsmsSpectrum]
        A dictionary with as keys the scan numbers and as values the
        corresponding spectra.
    """
    spectra = {}
    for spectrum_dict in mgf.read(mgf_filename):
        # TODO: Make sure the USI doesn't contain a peptide identification.
        cluster, usi = spectrum_dict['params']['title'].split(';')
        spectrum = sus.MsmsSpectrum(
            usi,
            spectrum_dict['params']['pepmass'][0],
            spectrum_dict['params']['charge'][0],
            spectrum_dict['m/z array'],
            spectrum_dict['intensity array'],
            retention_time=spectrum_dict['params']['rtinseconds'])
        spectrum.cluster = cluster
        if usi in spectra:
            raise ValueError(f'Non-unique USI: {usi}')
        spectra[usi] = spectrum
    return spectra


def get_scores(score_filename: str) -> pd.Series:
    """
    Read the PSM scores from the given MaxQuant identifications file.

    Parameters
    ----------
    score_filename : str
        The file name of the MaxQuant identifications file (msms.txt).

    Returns
    -------
    pd.Series
        A `Series` with as index the USI and as values the MaxQuant
        identification scores.
    """
    scores = pd.read_csv(score_filename, sep='\t',
                         usecols=['Raw file', 'Scan number', 'Score'])
    # FIXME: Don't hardcode the PXD identifier.
    scores['usi'] = ('mzspec:PXD004732:' + scores['Raw file'] +
                     '.raw::scan:' + scores['Scan number'].astype(str))
    scores = scores.set_index('usi')
    return scores['Score'].sort_index()


def get_best_representative(spectra: Dict[str, sus.MsmsSpectrum],
                            scores: pd.Series) -> sus.MsmsSpectrum:
    """
    Get the best representative spectrum with the highest score among the given
    spectra.

    If none of the given spectra have a score a ValueError is raised.

    If two spectra have the same maximum score, the spectrum with the
    alphanumeric first file name with the lowest scan number is selected.

    Parameters
    ----------
    spectra : Dict[str, sus.MsmsSpectrum]
        A dictionary with as keys the USIs and as values the corresponding
        spectra in a cluster.
    scores : pd.Series
        A `Series` with as index the USI and as values the identification
        scores.

    Returns
    -------
    sus.MsmsSpectrum
        The representative spectrum with the highest score.

    Raises
    ------
    ValueError
        If none of the given spectra have a score.
    """
    scores = scores[scores.index.isin(spectra.keys())]
    if len(scores) == 0:
        raise ValueError('No scores found for the given scan numbers')
    return spectra[scores.idxmax()]


def write_mgf(filename: str, spectra: List[sus.MsmsSpectrum]) -> None:
    """
    Write the given spectra to an MGF file.

    Parameters
    ----------
    filename : str
        The file name of the MGF output file.
    spectra : List[sus.MsmsSpectrum]
        The spectra to be written to the MGF file.
    """
    spectra_dict = [{'m/z array': spectrum.mz,
                     'intensity array': spectrum.intensity,
                     'params': {'title': (f'{spectrum.cluster};'
                                          f'{spectrum.identifier}'),
                                'pepmass': spectrum.precursor_mz,
                                'rtinseconds': spectrum.retention_time,
                                'charge': spectrum.precursor_charge}}
                    for spectrum in spectra]
    with open(filename, 'w') as f_out:
        mgf.write(spectra_dict, f_out)


def split_into_clusters(spectra: Dict[str, sus.MsmsSpectrum]) \
        -> Iterable[Dict[str, sus.MsmsSpectrum]]:
    """
    Yield collections of spectra grouped by cluster membership.

    Parameters
    ----------
    spectra : Dict[str, sus.MsmsSpectrum]
        The spectra to be grouped by cluster membership. The `MsmsSpectrum`
        objects should have a member variable `cluster` to designate which
        cluster they belong to.

    Returns
    -------
    Iterable[Dict[str, sus.MsmsSpectrum]]
        A dictionary with as keys the USI and as values the spectra for each
        individual cluster.
    """
    clusters = collections.defaultdict(list)
    for spectrum in spectra.values():
        clusters[spectrum.cluster].append(spectrum.identifier)
    for cluster_members in clusters.values():
        yield {usi: spectra[usi] for usi in cluster_members}


def best_spectrum(mgf_in_filename: str, mgf_out_filename: str,
                  scores_filename: str) -> None:
    """
    Represent clusters by their highest scoring member.

    Parameters
    ----------
    mgf_in_filename : str
        The file name of the MGF input file that contains the original spectra.
    mgf_out_filename : str
        The file name of the MGF output file that contains the cluster
        representatives.
    scores_filename : str
        The file name of the MaxQuant identifications file containing PSM
        scores.
    """
    scores = get_scores(scores_filename)
    spectra = get_cluster_spectra(mgf_in_filename)
    representatives = []
    for cluster in split_into_clusters(spectra):
        try:
            representatives.append(get_best_representative(cluster, scores))
        except ValueError:
            pass
    write_mgf(mgf_out_filename, representatives)


if __name__ == '__main__':
    best_spectrum(sys.argv[1], sys.argv[2], sys.argv[3])
