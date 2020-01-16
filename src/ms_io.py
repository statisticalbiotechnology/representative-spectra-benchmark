import os
from typing import Dict, Iterable

import pandas as pd
import pyteomics.mgf
import spectrum_utils.spectrum as sus
import tqdm


def read_spectra(filename: str) -> Iterable[sus.MsmsSpectrum]:
    ext = os.path.splitext(filename.lower())[1]
    if ext == '.mgf':
        yield from _read_spectra_mgf(filename)
    elif ext == '.mzml':
        yield from _read_spectra_mzml(filename)
    elif ext == '.mzxml':
        yield from _read_spectra_mzxml(filename)
    else:
        raise ValueError('Unsupported peak file format (supported formats: '
                         'MGF, mzML, mzXML)')


def _read_spectra_mgf(filename: str) -> Iterable[sus.MsmsSpectrum]:
    for spectrum_dict in tqdm.tqdm(pyteomics.mgf.read(filename),
                                   desc='Spectra read', unit='spectra'):
        yield sus.MsmsSpectrum(
            spectrum_dict['params']['title'],
            spectrum_dict['params']['pepmass'][0],
            spectrum_dict['params']['charge'][0],
            spectrum_dict['m/z array'],
            spectrum_dict['intensity array'],
            None,
            spectrum_dict['params']['rtinseconds'])


def _read_spectra_mzml(filename: str) -> Iterable[sus.MsmsSpectrum]:
    raise NotImplementedError


def _read_spectra_mzxml(filename: str) -> Iterable[sus.MsmsSpectrum]:
    raise NotImplementedError


###############################################################################


def read_clusters(filename: str, fmt: str) -> Dict[str, int]:
    if fmt == 'maracluster':
        return _read_clusters_maracluster(filename)
    elif fmt == 'pride-cluster':
        return _read_clusters_pridecluster(filename)
    elif fmt == 'ms-cluster':
        return _read_clusters_mscluster(filename)
    else:
        raise ValueError('Unsupported cluster file format (supported formats: '
                         'MaRaCluster, pride-cluster, MS-Cluster)')


def _read_clusters_maracluster(filename: str) -> Dict[str, int]:
    with open(filename) as f_in:
        clusters = {}
        cluster_i = 0
        for line in f_in:
            if not line.strip():
                cluster_i += 1
            else:
                filename, scan, *_ = line.split('\t')
                clusters[f'{os.path.splitext(filename)[0]}:scan:{scan}'] = \
                    cluster_i
        return clusters


def _read_clusters_pridecluster(filename: str) -> Dict[str, int]:
    raise NotImplementedError


def _read_clusters_mscluster(filename: str) -> Dict[str, int]:
    raise NotImplementedError


###############################################################################


def write_spectra(filename: str, spectra: Iterable[sus.MsmsSpectrum]) -> None:
    """
    Write the given spectra to an MGF file.

    Parameters
    ----------
    filename : str
        The file name of the MGF output file.
    spectra : List[Dict]
        The spectra as Pyteomics dictionaries to be written to the MGF file.
    """
    def _spectra_to_dicts(spectra: Iterable[sus.MsmsSpectrum]) \
            -> Iterable[Dict]:
        for spectrum in tqdm.tqdm(spectra, desc='Spectra written',
                                  unit='spectra'):
            yield {'params': {'title': spectrum.identifier,
                              'pepmass': spectrum.precursor_mz,
                              'rtinseconds': spectrum.retention_time,
                              'charge': spectrum.precursor_charge},
                   'm/z array': spectrum.mz,
                   'intensity array': spectrum.intensity}

    with open(filename, 'w') as f_out:
        pyteomics.mgf.write(_spectra_to_dicts(spectra), f_out)


###############################################################################


def read_maxquant_psms(msms_filename: str, px_accession: str) -> pd.Series:
    """
    Read the PSM scores from the given MaxQuant identifications file.

    Parameters
    ----------
    msms_filename : str
        The file name of the MaxQuant identifications file (msms.txt).
    px_accession : str
        PXD identifier for the dataset the identifications correspond to.
        Needed to compile the USI.

    Returns
    -------
    pd.Series
        A `DataFrame` with as index the USI and as values the MaxQuant PSM
        scores ("Score") and unmodified peptide sequences ("Sequence").
    """
    # TODO: Deal with the case when there are multiple PSMs per scan.
    scores = pd.read_csv(
        msms_filename, sep='\t', usecols=['Raw file', 'Scan number',
                                          'Sequence', 'Score'])
    scores['usi'] = (f'mzspec:{px_accession}:' + scores['Raw file'] +
                     ':scan:' + scores['Scan number'].astype(str))
    scores = scores[['usi', 'Sequence', 'Score']].set_index('usi')
    return scores.sort_index()
