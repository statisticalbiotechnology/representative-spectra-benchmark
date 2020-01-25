import json
import logging
import os
from typing import Dict, Iterable

import pandas as pd
import pyteomics.mgf
import pyteomics.mzml
import pyteomics.mzxml
import spectrum_utils.spectrum as sus
import tqdm
from lxml.etree import LxmlError


logger = logging.getLogger('cluster_representative')


def read_spectra(filename: str) -> Iterable[sus.MsmsSpectrum]:
    """
    Read MS/MS spectra from a peak file.

    Supported formats: MGF, mzML, mzXML.

    Parameters
    ----------
    filename : str
        The peak file name.

    Returns
    -------
    Iterable[sus.MsmsSpectrum]
        An iterable of spectra in the given peak file.
    """
    ext = os.path.splitext(filename.lower())[1]
    if ext == '.mgf':
        yield from _read_spectra_mgf(filename)
    elif ext == '.mzml':
        yield from _read_spectra_mzml(filename)
    elif ext == '.mzxml':
        yield from _read_spectra_mzxml(filename)
    else:
        logger.error('Unsupported peak file format (supported formats: MGF, '
                     'mzML, mzXML)')
        raise ValueError('Unsupported peak file format (supported formats: '
                         'MGF, mzML, mzXML)')


def _read_spectra_mgf(filename: str) -> Iterable[sus.MsmsSpectrum]:
    """
    Read MS/MS spectra from an MGF file.

    Parameters
    ----------
    filename : str
        The MGF file name.

    Returns
    -------
    Iterable[sus.MsmsSpectrum]
        An iterable of spectra in the given MGF file.
    """
    for spectrum_dict in tqdm.tqdm(pyteomics.mgf.read(filename),
                                   desc='Spectra read', unit='spectra'):
        spectrum = sus.MsmsSpectrum(
            spectrum_dict['params']['title'],
            spectrum_dict['params']['pepmass'][0],
            spectrum_dict['params']['charge'][0],
            spectrum_dict['m/z array'],
            spectrum_dict['intensity array'],
            None,
            spectrum_dict['params'].get('rtinseconds'))
        spectrum.filename = spectrum_dict['params'].get(
            'filename', os.path.splitext(os.path.basename(filename))[0])
        if 'scans' in spectrum_dict['params']:
            spectrum.scan = spectrum_dict['params']['scans']
        if 'cluster' in spectrum_dict['params']:
            spectrum.cluster = spectrum_dict['params']['cluster']
        yield spectrum


def _read_spectra_mzml(filename: str) -> Iterable[sus.MsmsSpectrum]:
    """
    Read MS/MS spectra from an mzML file.

    Parameters
    ----------
    filename : str
        The mzML file name.

    Returns
    -------
    Iterable[sus.MsmsSpectrum]
        An iterable of spectra in the given mzML file.
    """
    with pyteomics.mzml.MzML(filename) as f_in:
        try:
            for spectrum_dict in f_in:
                if int(spectrum_dict.get('ms level', -1)) == 2:
                    precursor = spectrum_dict['precursorList']['precursor'][0]
                    precursor_ion = (precursor['selectedIonList']
                                     ['selectedIon'][0])
                    if 'charge state' in precursor_ion:
                        precursor_charge = int(precursor_ion['charge state'])
                    elif 'possible charge state' in precursor_ion:
                        precursor_charge = int(
                            precursor_ion['possible charge state'])
                    else:
                        logger.warning('Unknown precursor charge, skipped '
                                       'spectrum...')
                        continue
                    spectrum = sus.MsmsSpectrum(
                        spectrum_dict['id'],
                        precursor_ion['selected ion m/z'],
                        precursor_charge,
                        spectrum_dict['m/z array'],
                        spectrum_dict['intensity array'],
                        None,
                        (spectrum_dict['scanList']['scan'][0]
                                      ['scan start time']))
                    if 'scan=' in spectrum.identifier:
                        spectrum.scan = int(
                            spectrum.identifier[
                                spectrum.identifier.find('scan=')
                                + len('scan='):])
                        spectrum.filename = os.path.splitext(
                            os.path.basename(filename))[0]
                    yield spectrum
        except LxmlError as e:
            logger.error('Failed to read file %s: %s', filename, e)


def _read_spectra_mzxml(filename: str) -> Iterable[sus.MsmsSpectrum]:
    """
    Read MS/MS spectra from an mzXML file.

    Parameters
    ----------
    filename : str
        The mzXML file name.

    Returns
    -------
    Iterable[sus.MsmsSpectrum]
        An iterable of spectra in the given mzXML file.
    """
    with pyteomics.mzxml.MzXML(filename) as f_in:
        try:
            for spectrum_dict in f_in:
                if int(spectrum_dict.get('msLevel', -1)) == 2:
                    if 'precursorCharge' in spectrum_dict['precursorMz'][0]:
                        precursor_charge = (spectrum_dict['precursorMz'][0]
                                                         ['precursorCharge'])
                    else:
                        logger.warning('Unknown precursor charge, skipped '
                                       'spectrum...')
                        continue
                    spectrum = sus.MsmsSpectrum(
                        spectrum_dict['id'],
                        spectrum_dict['precursorMz'][0]['precursorMz'],
                        precursor_charge,
                        spectrum_dict['m/z array'],
                        spectrum_dict['intensity array'],
                        None,
                        spectrum_dict['retentionTime'])
                    spectrum.scan = int(spectrum_dict['id'])
                    spectrum.filename = os.path.splitext(
                        os.path.basename(filename))[0]
                    yield spectrum
        except LxmlError as e:
            logger.warning('Failed to read file %s: %s', filename, e)


###############################################################################


def read_clusters(filename: str, fmt: str) -> Dict[str, int]:
    """
    Read cluster assignments.

    Parameters
    ----------
    filename : str
        The cluster assignment file name.
    fmt : str
        Format of the cluster assignment file. Supported formats:
        "maracluster", "spectra-cluster", "ms-cluster".

    Returns
    -------
    Dict[str, int]
        A dictionary with as keys the spectrum identifiers (format
        "{filename}:scan:{scan}") and as value the cluster index.
    """
    if fmt == 'maracluster':
        return _read_clusters_maracluster(filename)
    elif fmt == 'spectra-cluster':
        return _read_clusters_spectracluster(filename)
    elif fmt == 'ms-cluster':
        return _read_clusters_mscluster(filename)
    else:
        raise ValueError('Unsupported cluster file format (supported formats: '
                         'MaRaCluster, spectra-cluster, MS-Cluster)')


def _read_clusters_maracluster(filename: str) -> Dict[str, int]:
    """
    Read MaRaCluster cluster assignments.

    Parameters
    ----------
    filename : str
        The MaRaCluster cluster assignment file name.

    Returns
    -------
    Dict[str, int]
        A dictionary with as keys the spectrum identifiers (format
        "{filename}:scan:{scan}") and as value the cluster index.
    """
    with open(filename) as f_in:
        clusters, cluster_i = {}, 0
        for line in f_in:
            if not line.strip():
                cluster_i += 1
            else:
                fn, scan, *_ = line.split('\t')
                clusters[f'{os.path.splitext(fn)[0]}:scan:{scan}'] = \
                    cluster_i
        return clusters


def _read_clusters_spectracluster(filename: str) -> Dict[str, int]:
    """
    Read spectra-cluster cluster assignments.

    Parameters
    ----------
    filename : str
        The spectra-cluster cluster assignment file name.

    Returns
    -------
    Dict[str, int]
        A dictionary with as keys the spectrum identifiers (format
        "{filename}:scan:{scan}") and as value the cluster index.
    """
    with open(filename) as f_in:
        clusters, cluster_i = {}, 0
        for line in f_in:
            if line.startswith('=Cluster='):
                cluster_i += 1
            elif line.startswith('SPEC'):
                fn = line[line.find('#file') + len('#file'):line.find('#id')]
                scan = line[line.find('#id=index=') + len('#id=index='):
                            line.find('#title')]
                clusters[f'{os.path.splitext(fn)[0]}:scan:{scan}'] = \
                    cluster_i
        return clusters


def _read_clusters_mscluster(filename: str) -> Dict[str, int]:
    """
    Read MS-Cluster cluster assignments.

    Parameters
    ----------
    filename : str
        The MS-Cluster cluster assignment file name.

    Returns
    -------
    Dict[str, int]
        A dictionary with as keys the spectrum identifiers (format
        "{filename}:scan:{scan}") and as value the cluster index.
    """
    logger.warning('MS-Cluster output reading is not fully supported')
    with open(filename) as f_in:
        clusters, cluster_i = {}, 0
        for line in f_in:
            if line.startswith('mscluster'):
                cluster_i += 1
            elif not line.isspace():
                # FIXME: Missing file information for MS-Cluster output?
                clusters[int(line.split('\t')[2])] = cluster_i
        return clusters


###############################################################################


def write_spectra(filename: str, spectra: Iterable[sus.MsmsSpectrum]) -> None:
    def _spectra_to_dicts(spectra: Iterable[sus.MsmsSpectrum]) \
            -> Iterable[Dict]:
        for spectrum in tqdm.tqdm(spectra, desc='Spectra written',
                                  unit='spectra'):
            params = {'title': spectrum.identifier,
                      'pepmass': spectrum.precursor_mz,
                      'charge': spectrum.precursor_charge}
            if hasattr(spectrum, 'retention_time'):
                params['rtinseconds'] = spectrum.retention_time
            if hasattr(spectrum, 'filename'):
                params['filename'] = spectrum.filename
            if hasattr(spectrum, 'scan'):
                params['scans'] = spectrum.scan
            if hasattr(spectrum, 'cluster'):
                params['cluster'] = spectrum.cluster
            yield {'params': params,
                   'm/z array': spectrum.mz,
                   'intensity array': spectrum.intensity}

    with open(filename, 'w') as f_out:
        pyteomics.mgf.write(_spectra_to_dicts(spectra), f_out)


###############################################################################


def read_psms(filename: str) -> pd.DataFrame:
    ext = os.path.splitext(filename.lower())[1]
    if ext == '.mztab':
        return _read_psms_mztab(filename)
    elif ext == '.mzidentml':
        return _read_psms_mzidentml(filename)
    elif ext == '.json':
        return _read_psms_json(filename)
    elif os.path.basename(filename) == 'msms.txt':
        return _read_psms_maxquant(filename)
    else:
        raise ValueError('Unsupported PSM file format (supported formats: '
                         'mzTab, mzIdentML, JSON, MaxQuant)')


def _read_psms_mztab(filename: str) -> pd.DataFrame:
    raise NotImplementedError


def _read_psms_mzidentml(filename: str) -> pd.DataFrame:
    raise NotImplementedError


def _read_psms_json(filename: str) -> pd.DataFrame:
    raise NotImplementedError


def _read_psms_maxquant(filename: str) -> pd.DataFrame:
    psms = (pd.read_csv(filename, sep='\t', usecols=['Raw file', 'Scan number',
                                                     'Sequence', 'Score'])
            .rename(columns={'Sequence': 'sequence', 'Score': 'score'}))
    psms['spectra_ref'] = (psms['Raw file'] + ':scan:' +
                           psms['Scan number'].astype(str))
    return (psms[['spectra_ref', 'sequence', 'score']]
            .drop_duplicates('spectra_ref')
            .set_index('spectra_ref')
            .sort_index())


###############################################################################


def write_distance_dict_to_json(filename: str, distances: Dict) -> None:
    """
    Write the given distance dict to a JSON file.

    Parameters
    ----------
    filename : str
        The file name of the JSON output file.
    distances : Dict
        A dictionary where the keys are clusters and the values the distances.
    """
    # os.makedirs(RESULTFOLDER, exists_ok=True)
    # path = os.path.join(RESULTFOLDER, filename)
    with open(filename, 'w') as outfile:
        json.dump(distances, outfile, indent=4, sort_keys=True)
