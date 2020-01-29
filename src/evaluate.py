import functools
import logging

import click
import pandas as pd
import tqdm

import metrics
import ms_io
import representative


logger = logging.getLogger('cluster_representative')


@click.command('evaluate',
               help='Evaluate how well representative spectra represent '
                    'cluster members')
@click.option('--filename_spectra', 'filename_spectra',
              type=click.Path(dir_okay=False),
              help='Input MGF file containing original spectra (including '
                   'cluster assignments)')
@click.option('--filename_representatives', 'filename_representatives',
              type=click.Path(dir_okay=False),
              help='Input MGF file containing cluster representatives')
@click.option('--filename_out', 'filename_out',
              type=click.Path(dir_okay=False),
              help='Output CSV file containing representative scores')
@click.option('--measure', 'measure',
              type=click.Choice(['avg_dot', 'fraction_by']),
              help='Measure used to evaluate cluster representatives (options:'
                   ' "avg_dot", "fraction_by")')
@click.option('--fragment_mz_tolerance', 'fragment_mz_tolerance', type=float,
              default=0.02, show_default=True,
              help='Fragment m/z tolerance used during spectrum comparison or'
                   ' fragment annotation (optional; required for the "avg_dot"'
                   ' and "fraction_by" methods)')
@click.option('--filename_psm', 'filename_psm',
              type=click.Path(exists=True, dir_okay=False),
              help='Input PSM file (optional; supported formats: mzTab, '
                   'mzIdentML, JSON, MaxQuant; required for the '
                   '"fraction_by" method)')
def evaluate(filename_spectra: str, filename_representatives: str,
             filename_out: str, metric: str,
             fragment_mz_tolerance: float = 0.02,
             filename_psm: str = None) -> None:
    """
    Evaluate how well representative spectra represent cluster members.

    Parameters
    ----------
    filename_spectra : str
        Input MGF file containing original spectra (including cluster
        assignments).
    filename_representatives : str
        Input MGF file containing cluster representatives.
    filename_out : str
        Output JSON file containing representative scores.
    metric : str
        Measure used to evaluate cluster representatives (options: "avg_dot").
    fragment_mz_tolerance : float
        Fragment m/z tolerance used during spectrum comparison (optional;
        required for the "avg_dot" method).
    filename_psm : str
        Input PSM file (optional; supported formats: mzTab, mzIdentML, JSON,
        MaxQuant; required for the "fraction_by" method).
    """
    if metric == 'avg_dot':
        metric_func = functools.partial(
            metrics.avg_dot, fragment_mz_tolerance=fragment_mz_tolerance)
    elif metric == 'fraction_by':
        metric_func = functools.partial(
            metrics.fraction_by, fragment_mz_tolerance=fragment_mz_tolerance)
    else:
        raise ValueError('Unknown metric specified: %s', metric)

    logger.info('Read original spectra from spectrum file %s',
                filename_spectra)
    spectra = {f'{spectrum.filename}:scan:{spectrum.scan}': spectrum
               for spectrum in ms_io.read_spectra(filename_spectra)}
    if metric == 'fraction_by':
        logger.info('Assign peptide identifications from file %s to the '
                    'cluster representatives', filename_psm)
        psms = ms_io.read_psms(filename_psm)
        for spectra_ref, sequence in psms['sequence'].items():
            if spectra_ref in spectra:
                spectra[spectra_ref].peptide = sequence

    logger.info('Read cluster representatives from spectrum file %s',
                filename_representatives)
    representatives = {
        spectrum.cluster: spectrum
        for spectrum in ms_io.read_spectra(filename_representatives)}

    logger.info('Evaluate cluster representatives')
    cluster_keys, scores = [], []
    for cluster_key, cluster_spectra in tqdm.tqdm(
            representative.get_cluster_spectra(spectra),
            desc='Clusters processed', unit='clusters'):
        if cluster_key in representatives:
            cluster_keys.append(cluster_key)
            scores.append(metric_func(
                representatives[cluster_key], cluster_spectra.values()))

    logger.info('Export representative evaluation scores to JSON file %s',
                filename_out)
    (pd.DataFrame({'cluster': cluster_keys, 'score': scores})
     .to_csv(filename_out, index=False))


if __name__ == '__main__':
    logging.basicConfig(format='{asctime} [{levelname}/{processName}] '
                               '{module}.{funcName} : {message}',
                        style='{', level=logging.INFO)

    evaluate()

    logging.shutdown()
