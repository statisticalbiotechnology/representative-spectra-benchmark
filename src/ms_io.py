import os
import re
from typing import Dict, Iterable

import pandas as pd
import pyteomics.mgf
import spectrum_utils.spectrum as sus


def read_cluster_spectra(mgf_filename: str, usi_present: bool=True, cluster_present:bool=True) -> Dict[str, sus.MsmsSpectrum]:
    """
    Read all spectra with cluster information from the given MGF file.

    Parameters
    ----------
    mgf_filename : str
        The file name of the MGF file to be read.
    usi_present: bool
        Do the titles of the mgf contain a USI
    cluster_present: bool
        Do the titles of the mgf contain a cluster ID

    Raises
    -------
    ValueError:
        In case of duplicate Identifiers (being USI or Cluster IDs) a 
        ValueError is raised
    NotImplementedError: 
        In case only USIs are present, the functionality is yet undefined.


    Returns
    -------
    Dict[str, sus.MsmsSpectrum]
        A dictionary with keys the USI (without optional peptide
        identification) and as values the corresponding spectra
        or a dictionary with clusters as keys in case no USI is 
        given in `mgf_filename`.
    """
    spectra = {}
    for spectrum_dict in pyteomics.mgf.read(mgf_filename):
        spectrum = _dict_to_spectrum(spectrum_dict)
        if usi_present and cluster_present:
            title = re.match(r'(cluster-\d+);(mzspec:\w+:.+:(scan|index):\d+)',
                         spectrum.identifier)
            spectrum.cluster, spectrum.identifier = title.group(1), title.group(2)
            if spectrum.identifier in spectra:
                raise ValueError(f'Non-unique USI: {spectrum.identifier}')
            spectra[spectrum.identifier] = spectrum
        elif cluster_present and not usi_present:
            title = re.match(r'(cluster-\d+)',
                         spectrum.identifier)
            spectrum.cluster = title.group(1)    
            spectra[spectrum.cluster] = spectrum 
            if spectrum.cluster in spectra:
                raise ValueError(f'Non unique cluster identifier: {spectrum.cluster}')
        else:
            raise NotImplementedError("Missing functionality for now.")
    return spectra


def _dict_to_spectrum(spectrum_dict: Dict):
    return sus.MsmsSpectrum(
        spectrum_dict['params']['title'],
        spectrum_dict['params']['pepmass'][0],
        spectrum_dict['params']['charge'][0],
        spectrum_dict['m/z array'],
        spectrum_dict['intensity array'],
        retention_time=spectrum_dict['params']['rtinseconds'])


def _spectrum_to_dict(spectrum: sus.MsmsSpectrum):
    return {'m/z array': spectrum.mz,
            'intensity array': spectrum.intensity,
            'params': {'title': spectrum.identifier,
                       'pepmass': spectrum.precursor_mz,
                       'rtinseconds': spectrum.retention_time,
                       'charge': spectrum.precursor_charge}}


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


def read_maracluster_clusters(maracluster_filename: str, px_accession: str):
    with open(maracluster_filename) as f_in:
        usis, clusters = [], []
        cluster_i = 0
        for line in f_in:
            if not line.strip():
                cluster_i += 1
            else:
                filename, scan, *_ = line.split('\t')
                usis.append(_build_usi(
                    px_accession, os.path.splitext(filename)[0], int(scan)))
                clusters.append(cluster_i)
        return pd.Series(clusters, usis, name='cluster')


def _build_usi(px_accession: str, raw_name: str, scan: int,
               peptide_sequence: str = None, charge: int = None,
               cluster_id: str = None):
    usi = f'mzspec:{px_accession}:{raw_name}:scan:{scan}'
    if peptide_sequence is not None and charge is not None:
        usi = f'{usi}:{peptide_sequence}/{charge}'
    if cluster_id is not None:
        usi = f'{cluster_id};{usi}'
    return usi


def write_mgf(filename: str, spectra: Iterable[Dict]) -> None:
    """
    Write the given spectra to an MGF file.

    Parameters
    ----------
    filename : str
        The file name of the MGF output file.
    spectra : List[Dict]
        The spectra as Pyteomics dictionaries to be written to the MGF file.
    """
    with open(filename, 'w') as f_out:
        pyteomics.mgf.write(spectra, f_out)
