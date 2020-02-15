import collections
import logging
from typing import Dict, Iterable

import click
import spectrum_utils.spectrum as sus
import tqdm

import ms_io
import selector


logger = logging.getLogger('cluster_representative')


def get_cluster_spectra(spectra: Dict[str, sus.MsmsSpectrum]) \
        -> Iterable[Dict[str, sus.MsmsSpectrum]]:
    """
    Get collections of spectra for each cluster.

    Parameters
    ----------
    spectra : Dict[str, sus.MsmsSpectrum]
        A dictionary of spectra with as keys the spectrum identifiers and as
        values the MsmsSpectrum object.

    Returns
    -------
    Iterable[Dict[str, sus.MsmsSpectrum]]
        An iterable of spectrum dictionaries (see above) grouped by cluster.
    """
    clusters = collections.defaultdict(list)
    for spectrum_key, spectrum in spectra.items():
        clusters[spectrum.cluster].append(spectrum_key)
    for cluster_key, cluster_members in clusters.items():
        yield cluster_key, {spectrum_key: spectra[spectrum_key]
                            for spectrum_key in cluster_members}


@click.command('representative',
               help='Export representative spectra for each cluster')
@click.option('--filename_in', 'filename_in',
              type=click.Path(exists=True, dir_okay=False),
              help='Input spectrum file containing cluster assignments')
@click.option('--filename_out', 'filename_out',
              type=click.Path(dir_okay=False),
              help='Output spectrum file containing representative spectra for'
                   ' each cluster')
@click.option('--representative_method', 'representative_method',
              type=click.Choice(['most_similar', 'best_spectrum', 'bin']),
              help='Method used to select the representative spectrum for each'
                   ' cluster (options: "best_spectrum", "most_similar", '
                   '"bin")')
@click.option('--min_cluster_size', 'min_cluster_size', type=int,
              default=5, show_default=True,
              help='Only consider clusters consisting of a minimum number of '
                   'spectra (optional)')
# Options for BestSpectrumRepresentativeSelector.
@click.option('--filename_psm', 'filename_psm',
              type=click.Path(exists=True, dir_okay=False),
              help='Input PSM file (optional; supported formats: mzTab, '
                   'mzIdentML, JSON, MaxQuant; required for the '
                   '"best_spectrum" method)')
@click.option('--higher_is_better/--lower_is_better', 'higher_is_better',
              default=True, show_default=True,
              help='Flag indicating whether higher PSM scores are better or '
                   'not')
# Options for MostSimilarRepresentativeSelector.
@click.option('--sim', 'sim', type=click.Choice(['dot']),
              default='dot', show_default=True,
              help='Similarity measure to compare spectra to each other '
                   '(optional; required for the "most_similar" method; '
                   'options: "dot")')
@click.option('--fragment_mz_tolerance', 'fragment_mz_tolerance', type=float,
              default=0.02, show_default=True,
              help='Fragment m/z tolerance used during spectrum comparison '
                   '(optional; required for the "most_similar" method)')
# Options for BinningRepresentativeSelector.
@click.option('--min_mz', 'min_mz', type=float,
              default=100., show_default=True,
              help='Minimum m/z to consider for spectrum binning (optional; '
                   'required for the "bin" method)')
@click.option('--max_mz', 'max_mz', type=float,
              default=2000., show_default=True,
              help='Maximum m/z to consider for spectrum binning (optional; '
                   'required for the "bin" method)')
@click.option('--bin_size', 'bin_size', type=float,
              default=0.02, show_default=True,
              help='Bin size in m/z used for spectrum binning (optional; '
                   'required for the "bin" method)')
@click.option('--bin_peak_quorum', 'bin_peak_quorum', type=float,
              default=0.25, show_default=True,
              help='Relative number of spectra in a cluster that need to '
                   'contain a peak for it to be included in the representative'
                   ' spectrum (optional; required for the "bin" method)')
@click.option('--bin_edge_case_threshold', 'bin_edge_case_threshold',
              type=float, default=0.5, show_default=True,
              help='During binning try to correct m/z edge cases where the '
                   'm/z is closer to the bin edge than the given relative bin '
                   'size threshold (optional; required for the "bin" method)')
def representative(filename_in: str, filename_out: str,
                   representative_method: str, min_cluster_size: int = 5,
                   filename_psm: str = None, higher_is_better: bool = True,
                   sim: str = 'dot', fragment_mz_tolerance: float = 0.02,
                   min_mz: float = 100., max_mz: float = 2000.,
                   bin_size: float = 0.02, bin_peak_quorum: float = 0.25,
                   bin_edge_case_threshold: float = 0.5) -> None:
    """
    Select representative spectra for clusters using a specific selection
    method.

    Parameters
    ----------
    filename_in : str
        Input spectrum file containing cluster assignments.
    filename_out : str
        Output spectrum file containing representative spectra for each cluster.
    representative_method : str
        Method used to select the representative spectrum for each cluster
        (options: "best_spectrum", "most_similar", "bin").
    min_cluster_size : int
        Only consider clusters consisting of a minimum number of spectra
        (optional).
    filename_psm : str
        Input PSM file (optional; supported formats: mzTab, mzIdentML, JSON,
        MaxQuant; required for the "best_spectrum" method).
    higher_is_better : bool
        Flag indicating whether higher PSM scores are better or not.
    sim : str
        Similarity measure to compare spectra to each other (optional; required
        for the "most_similar" method; options: "dot").
    fragment_mz_tolerance : float
        Fragment m/z tolerance used during spectrum comparison (optional;
        required for the "most_similar" method).
    min_mz : float
        Minimum m/z to consider for spectrum binning (optional; required for
        the "bin" method).
    max_mz : float
        Maximum m/z to consider for spectrum binning (optional; required for
        the "bin" method).
    bin_size : float
        Bin size in m/z used for spectrum binning (optional; required for the
        "bin" method).
    bin_peak_quorum : float
        Relative number of spectra in a cluster that need to contain a peak for
        it to be included in the representative spectrum (optional; required
        for the "bin" method).
    bin_edge_case_threshold : float
        During binning try to correct m/z edge cases where the m/z is closer to
        the bin edge than the given relative bin size threshold (optional;
        required for the "bin" method).
    """
    logger.info('Read spectra from spectrum file %s', filename_in)
    spectra = {f'{spectrum.filename}:scan:{spectrum.scan}': spectrum
               for spectrum in ms_io.read_spectra(filename_in)}

    if representative_method == 'best_spectrum':
        rs = selector.BestSpectrumRepresentativeSelector(
            filename_psm, higher_is_better)
        logger.info('Select cluster representatives using the best spectrum '
                    'method')
    elif representative_method == 'most_similar':
        rs = selector.MostSimilarRepresentativeSelector(
            sim, fragment_mz_tolerance)
        logger.info('Select cluster representatives using the most similar '
                    'spectrum method')
    elif representative_method == 'bin':
        rs = selector.BinningRepresentativeSelector(
            min_mz, max_mz, bin_size, bin_peak_quorum, bin_edge_case_threshold)
        logger.info('Select cluster representatives via binning')
    else:
        raise ValueError('Unknown method to select the representative spectra')

    cluster_representatives = []
    for cluster_key, cluster_spectra in tqdm.tqdm(
            get_cluster_spectra(spectra), desc='Clusters processed',
            unit='clusters'):
        try:
            if len(cluster_spectra) >= min_cluster_size:
                cluster_representative = rs.select_representative(
                    cluster_spectra)
                cluster_representative.identifier = \
                    f'{rs.get_description()}:cluster={cluster_key}'
                cluster_representatives.append(cluster_representative)
        except ValueError:
            # logger.warning('No representative found for cluster %s',
            #                cluster_key)
            pass

    logger.info('Export %d cluster representatives to spectrum file %s',
                len(cluster_representatives), filename_out)
    ms_io.write_spectra(filename_out, cluster_representatives)


if __name__ == '__main__':
    logging.basicConfig(format='{asctime} [{levelname}/{processName}] '
                               '{module}.{funcName} : {message}',
                        style='{', level=logging.INFO)

    representative()

    logging.shutdown()
