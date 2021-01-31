import logging
import os
from typing import Dict, Iterable, List

import pandas as pd
import pyopenms
import pyteomics.mgf
import pyteomics.mzid
import pyteomics.mzml
import pyteomics.mztab
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
        if 'scan' in spectrum_dict['params']:
            spectrum.scan = spectrum_dict['params']['scan']
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
            for spectrum_dict in tqdm.tqdm(f_in, desc='Spectra read',
                                           unit='spectra'):
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
                    spectrum.filename = spectrum_dict.get(
                        'filename',
                        os.path.splitext(os.path.basename(filename))[0])
                    if 'scan' in spectrum_dict:
                        spectrum.scan = str(int(spectrum_dict['scan']))
                    elif 'scan=' in spectrum.identifier:
                        spectrum.scan = int(
                            spectrum.identifier[
                                spectrum.identifier.find('scan=')
                                + len('scan='):])
                    if 'cluster' in spectrum_dict:
                        spectrum.cluster = int(spectrum_dict['cluster'])
                    yield spectrum
        except LxmlError as e:
            logger.error('Failed to read file %s: %s', filename, e)


def _read_spectra_mzxml(filename: str) -> Iterable[sus.MsmsSpectrum]:
    """
    Read MS/MS spectra from an mzXML file.

    Attention: Reading intermediate mzXML files with clustering information is
    not supported, only original mzXML files.

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
            for spectrum_dict in tqdm.tqdm(f_in, desc='Spectra read',
                                           unit='spectra'):
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


def read_clusters(filename: str, fmt: str,
                  spectra_keys: List[str] = None) \
        -> Dict[str, int]:
    """
    Read cluster assignments.

    Parameters
    ----------
    filename : str
        The cluster assignment file name.
    fmt : str
        Format of the cluster assignment file. Supported formats:
        "maracluster", "spectra-cluster", "ms-cluster".
    spectra_keys : List[str]
        Required for MS-Cluster cluster assignment reading, ignored otherwise.
        Ordered list of spectrum keys used to properly set the spectrum
        identifiers from the MS-Cluster output.

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
        if spectra_keys is None:
            raise ValueError('The corresponding spectrum keys need to be '
                             'provided for MS-Cluster cluster assignment '
                             'reading')
        return _convert_clusters_mscluster(
            _read_clusters_mscluster(filename), spectra_keys)
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
    with open(filename) as f_in, tqdm.tqdm(
            desc='Cluster assignments read', unit='spectra') as progress_bar:
        clusters, cluster_i = {}, 0
        for line in f_in:
            if not line.strip():
                cluster_i += 1
            else:
                fn, scan, *_ = line.split('\t')
                clusters[f'{os.path.splitext(fn)[0]}:scan:{scan}'] = \
                    cluster_i
                progress_bar.update(1)
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
    with open(filename) as f_in, tqdm.tqdm(
            desc='Cluster assignments read', unit='spectra') as progress_bar:
        clusters, cluster_i = {}, 0
        for line in f_in:
            if line.startswith('=Cluster='):
                cluster_i += 1
            elif line.startswith('SPEC'):
                fn = line[line.find('#file') + len('#file') + 1:
                          line.find('#id')]
                scan = line[line.find('#id=index=') + len('#id=index='):
                            line.find('#title')]
                clusters[f'{os.path.splitext(os.path.basename(fn))[0]}:scan:'
                         f'{scan}'] = cluster_i
                progress_bar.update(1)
        return clusters


def _read_clusters_mscluster(filename: str) -> Dict[int, int]:
    """
    Read MS-Cluster cluster assignments.

    Parameters
    ----------
    filename : str
        The MS-Cluster cluster assignment file name.

    Returns
    -------
    Dict[int, int]
        A dictionary with as keys the spectrum identifiers (format
        "{filename}:scan:{scan}") and as value the cluster index.
    """
    logger.warning('MS-Cluster output reading is not fully supported')
    with open(filename) as f_in, tqdm.tqdm(
            desc='Cluster assignments read', unit='spectra') as progress_bar:
        clusters, cluster_i = {}, 0
        for line in f_in:
            if line.startswith('mscluster'):
                cluster_i += 1
            elif not line.isspace():
                clusters[int(line.split('\t')[2])] = cluster_i
                progress_bar.update(1)
        return clusters


def _convert_clusters_mscluster(clusters: Dict[int, int],
                                spectra_keys: List[str]) \
        -> Dict[str, int]:
    """
    Associate spectrum identifiers with MS-Cluster cluster assignments.

    Parameters
    ----------
    clusters : Dict[int, int]
        A dictionary with as keys the spectrum index and as value the cluster
        index.
    spectra_keys : Iterable[str]
        Ordered list of spectrum keys (format "{filename}:scan:{scan}").

    Returns
    -------
    Dict[str, int]
        A dictionary with as keys the spectrum identifiers (format
        "{filename}:scan:{scan}") and as value the cluster index.
    """
    return {spectra_keys[key_index]: cluster_i
            for key_index, cluster_i in clusters.items()}


###############################################################################


def write_spectra(filename: str, spectra: Iterable[sus.MsmsSpectrum]) -> None:
    """
    Write the given spectra to a peak file.

    Supported formats: MGF, mzML.

    Parameters
    ----------
    filename : str
        The file name where the spectra will be written.
    spectra : Iterable[sus.MsmsSpectrum]
        The spectra to be written to the peak file.
    """
    ext = os.path.splitext(filename.lower())[1]
    if ext == '.mgf':
        _write_spectra_mgf(filename, spectra)
    elif ext == '.mzml':
        _write_spectra_mzml(filename, spectra)
    else:
        logger.error('Unsupported peak file format (supported formats: MGF, '
                     'mzML)')
        raise ValueError('Unsupported peak file format (supported formats: '
                         'MGF, mzML)')


def _write_spectra_mgf(filename: str, spectra: Iterable[sus.MsmsSpectrum]) \
        -> None:
    """
    Write the given spectra to an MGF file.

    Parameters
    ----------
    filename : str
        The MGF file name where the spectra will be written.
    spectra : Iterable[sus.MsmsSpectrum]
        The spectra to be written to the MGF file.
    """
    with open(filename, 'w') as f_out:
        pyteomics.mgf.write(_spectra_to_dicts(spectra), f_out)


def _spectra_to_dicts(spectra: Iterable[sus.MsmsSpectrum]) -> Iterable[Dict]:
    """
    Convert MsmsSpectrum objects to Pyteomics MGF spectrum dictionaries.

    Parameters
    ----------
    spectra : Iterable[sus.MsmsSpectrum]
        The spectra to be converted to Pyteomics MGF dictionaries.

    Returns
    -------
    Iterable[Dict]
        The given spectra as Pyteomics MGF dictionaries.
    """
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
            params['scan'] = spectrum.scan
        if hasattr(spectrum, 'cluster'):
            params['cluster'] = spectrum.cluster
        yield {'params': params,
               'm/z array': spectrum.mz,
               'intensity array': spectrum.intensity}


def _write_spectra_mzml(filename: str, spectra: Iterable[sus.MsmsSpectrum]) \
        -> None:
    """
    Write the given spectra to an mzML file.

    Parameters
    ----------
    filename : str
        The mzML file name where the spectra will be written.
    spectra : Iterable[sus.MsmsSpectrum]
        The spectra to be written to the mzML file.
    """
    experiment = pyopenms.MSExperiment()
    for spectrum in tqdm.tqdm(spectra, desc='Spectra written', unit='spectra'):
        mzml_spectrum = pyopenms.MSSpectrum()
        mzml_spectrum.setMSLevel(2)
        mzml_spectrum.setNativeID(spectrum.identifier)
        precursor = pyopenms.Precursor()
        precursor.setMZ(spectrum.precursor_mz)
        precursor.setCharge(spectrum.precursor_charge)
        mzml_spectrum.setPrecursors([precursor])
        mzml_spectrum.set_peaks([spectrum.mz, spectrum.intensity])
        if hasattr(spectrum, 'retention_time'):
            mzml_spectrum.setRT(spectrum.retention_time)
        if hasattr(spectrum, 'filename'):
            mzml_spectrum.setMetaValue(
                'filename', str.encode(spectrum.filename))
        if hasattr(spectrum, 'scan'):
            mzml_spectrum.setMetaValue('scan', str.encode(str(spectrum.scan)))
        if hasattr(spectrum, 'cluster'):
            mzml_spectrum.setMetaValue(
                'cluster', str.encode(str(spectrum.cluster)))
        experiment.addSpectrum(mzml_spectrum)
    pyopenms.MzMLFile().store(filename, experiment)


###############################################################################


def read_psms(filename: str) -> pd.DataFrame:
    """
    Read spectrum identifications.

    Parameters
    ----------
    filename : str
        The PSM file name. Supported formats: mzTab, mzIdentML, JSON, MaxQuant.

    Returns
    -------
    pd.DataFrame
        A Pandas DataFrame with as rows the PSMs and as columns "sequence" and
        "score", indexed by their spectrum reference in the form of
        {filename}:scan:{scan}.
    """
    ext = os.path.splitext(filename.lower())[1]
    if ext == '.mztab':
        return _read_psms_mztab(filename)
    elif ext == '.mzidentml' or ext == '.mzid':
        return _read_psms_mzidentml(filename)
    elif ext == '.idxml':
        return _read_psms_idxml(filename)
    elif ext == '.json':
        return _read_psms_json(filename)
    elif os.path.basename(filename) == 'msms.txt':
        return _read_psms_maxquant(filename)
    else:
        raise ValueError('Unsupported PSM file format (supported formats: '
                         'mzTab, mzIdentML, JSON, MaxQuant)')


def _read_psms_mztab(filename: str) -> pd.DataFrame:
    """
    Read mzTab spectrum identifications.

    To correctly have the PSM index the "spectra_ref" column of the mzTab file
    should use the scan number specification.

    Parameters
    ----------
    filename : str
        The mzTab file name.

    Returns
    -------
    pd.DataFrame
        A Pandas DataFrame with as rows the PSMs and as columns "sequence" and
        "score", indexed by their spectrum reference in the form of
        {filename}:scan:{scan}.
    """
    mztab = pyteomics.mztab.MzTab(filename)
    run_names = {key[key.find('ms_run['):key.find('-location')]:
                     os.path.splitext(os.path.basename(value))[0]
                 for key, value in mztab.metadata.items()
                 if key.startswith('ms_run') and key.endswith('location')}
    psms = (mztab.spectrum_match_table
            .rename(columns={'search_engine_score[1]': 'score'}))
    if not all(psms['spectra_ref'].str.contains('scan')):
        raise ValueError('Spectrum references with scan information required')
    run = pd.Series([spectra_ref[:spectra_ref.find(':')]
                     for spectra_ref in psms['spectra_ref']])
    scan = pd.Series([spectra_ref[spectra_ref.find('=') + 1:]
                      for spectra_ref in psms['spectra_ref']])
    psms['spectra_ref'] = (run.map(run_names) + ':scan:' + scan).values
    return (psms[['spectra_ref', 'sequence', 'score']]
            .drop_duplicates('spectra_ref')
            .set_index('spectra_ref'))


def _read_psms_mzidentml(filename: str) -> pd.DataFrame:
    """
    Read mzIdentML spectrum identifications.

    Parameters
    ----------
    filename : str
        The mzIdentML file name.

    Returns
    -------
    pd.DataFrame
        A Pandas DataFrame with as rows the PSMs and as columns "sequence" and
        "score", indexed by their spectrum reference in the form of
        {filename}:scan:{scan}.
    """
    raise NotImplementedError('mzIdentML support is incomplete')
    filenames, scan, sequences, score = [], [], [], []
    for psm_dict in tqdm.tqdm(pyteomics.mzid.read(filename),
                              desc='PSMs read', unit='PSMs'):
        filenames.append(os.path.splitext(psm_dict['name'])[0])
        scan.append(-1)     # FIXME
        sequences.append(
            psm_dict['SpectrumIdentificationItem'][0]['PeptideSequence'])
        score.append(-1)    # FIXME
    psms = pd.DataFrame({'filename': filenames, 'scan': scan,
                         'sequence': sequences, 'score': score})
    psms['spectra_ref'] = (psms['filename'] + ':scan:' +
                           psms['scan'].astype(str))
    return (psms[['spectra_ref', 'sequence', 'score']]
            .drop_duplicates('spectra_ref')
            .set_index('spectra_ref'))


def _read_psms_idxml(filename: str) -> pd.DataFrame:
    """
    Read idXML spectrum identifications.

    Parameters
    ----------
    filename : str
        The idXML file name.

    Returns
    -------
    pd.DataFrame
        A Pandas DataFrame with as rows the PSMs and as columns "sequence" and
        "score", indexed by their spectrum reference in the form of
        {filename}:scan:{scan}.
    """
    protein_ids, psms, scans, sequences, scores = [], [], [], [], []
    pyopenms.IdXMLFile().load(filename, protein_ids, psms)
    peak_filename = os.path.splitext(os.path.basename(
        protein_ids[0].getMetaValue('spectra_data')[0].decode()))[0]
    for psm in tqdm.tqdm(psms, desc='PSMs read', unit='PSMs'):
        spectrum_index = psm.getMetaValue('spectrum_reference').decode()
        scans.append(
            int(spectrum_index[spectrum_index.find('scan=') + len('scan='):]))
        sequences.append(psm.getHits()[0].getSequence().toString().decode())
        scores.append(psm.getHits()[0].getScore())
    psms = pd.DataFrame({'filename': peak_filename, 'scan': scans,
                         'sequence': sequences, 'score': scores})
    psms['spectra_ref'] = (psms['filename'] + ':scan:' +
                           psms['scan'].astype(str))
    return (psms[['spectra_ref', 'sequence', 'score']]
            .drop_duplicates('spectra_ref')
            .set_index('spectra_ref'))


def _read_psms_json(filename: str) -> pd.DataFrame:
    """
    Read JSON spectrum identifications.

    Parameters
    ----------
    filename : str
        The JSON file name.

    Returns
    -------
    pd.DataFrame
        A Pandas DataFrame with as rows the PSMs and as columns "sequence" and
        "score", indexed by their spectrum reference in the form of
        {filename}:scan:{scan}.
    """
    raise NotImplementedError


def _read_psms_maxquant(filename: str) -> pd.DataFrame:
    """
    Read MaxQuant spectrum identifications.

    Parameters
    ----------
    filename : str
        The MaxQuant PSM file name (msms.txt).

    Returns
    -------
    pd.DataFrame
        A Pandas DataFrame with as rows the PSMs and as columns "sequence" and
        "score", indexed by their spectrum reference in the form of
        {filename}:scan:{scan}.
    """
    psms = (pd.read_csv(filename, sep='\t', usecols=['Raw file', 'Scan number',
                                                     'Sequence', 'Score'])
            .rename(columns={'Sequence': 'sequence', 'Score': 'score'}))
    psms['spectra_ref'] = (psms['Raw file'] + ':scan:' +
                           psms['Scan number'].astype(str))
    return (psms[['spectra_ref', 'sequence', 'score']]
            .drop_duplicates('spectra_ref')
            .set_index('spectra_ref'))
