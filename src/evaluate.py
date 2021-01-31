import functools
import logging
from typing import Tuple

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
              help='Input spectrum file containing original spectra'
                   '(including cluster assignments)')
@click.option('--filename_representatives', 'filename_representatives',
              type=click.Path(dir_okay=False),
              help='Input spectrum file containing cluster representatives')
@click.option('--filename_out', 'filename_out',
              type=click.Path(dir_okay=False),
              help='Output CSV file containing representative scores')
@click.option('--measure', 'measures',
              type=click.Choice(['avg_dot', 'fraction_by']), multiple=True,
              help='Measure(s) used to evaluate cluster representatives '
                   '(options: "avg_dot", "fraction_by")')
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
             filename_out: str, measures: Tuple[str],
             fragment_mz_tolerance: float = 0.02,
             filename_psm: str = None) -> None:
    """
    Evaluate how well representative spectra represent cluster members.

    Parameters
    ----------
    filename_spectra : str
        Input spectrum file containing original spectra (including cluster
        assignments).
    filename_representatives : str
        Input spectrum file containing cluster representatives.
    filename_out : str
        Output JSON file containing representative scores.
    measures : Tuple[str]
        Measure(s) used to evaluate cluster representatives (options:
        '"avg_dot", "fraction_by")).
    fragment_mz_tolerance : float
        Fragment m/z tolerance used during spectrum comparison (optional;
        required for the "avg_dot" method).
    filename_psm : str
        Input PSM file (optional; supported formats: mzTab, mzIdentML, JSON,
        MaxQuant; required for the "fraction_by" method).
    """
    measure_funcs, scores = {}, {}
    for measure in measures:
        if measure == 'avg_dot':
            measure_funcs[measure] = functools.partial(
                metrics.avg_dot, fragment_mz_tolerance=fragment_mz_tolerance)
            scores[measure] = []
        elif measure == 'fraction_by':
            measure_funcs[measure] = functools.partial(
                metrics.fraction_by,
                fragment_mz_tolerance=fragment_mz_tolerance)
            scores[measure] = []
        else:
            raise ValueError('Unknown measure specified: %s', measure)

    logger.info('Read original spectra from spectrum file %s',
                filename_spectra)
    spectra = {f'{spectrum.filename}:scan:{spectrum.scan}': spectrum
               for spectrum in ms_io.read_spectra(filename_spectra)}

    logger.info('Read cluster representatives from spectrum file %s',
                filename_representatives)
    representatives = {
        spectrum.cluster: spectrum
        for spectrum in ms_io.read_spectra(filename_representatives)}
    if 'fraction_by' in measures:
        logger.info('Assign peptide identifications from file %s to the '
                    'cluster representatives', filename_psm)
        psms = ms_io.read_psms(filename_psm)
        for spectra_ref, sequence in tqdm.tqdm(
                psms['sequence'].items(), desc='PSMs assigned', unit='PSMs'):
            if spectra_ref in representatives:
                representatives[spectra_ref].peptide = sequence

    logger.info('Evaluate cluster representatives')
    cluster_keys = []
    for cluster_key, cluster_spectra in tqdm.tqdm(
            representative.get_cluster_spectra(spectra),
            desc='Clusters processed', unit='clusters'):
        if cluster_key in representatives:
            cluster_keys.append(cluster_key)
            for measure, measure_func in measure_funcs.items():
                scores[measure].append(measure_func(
                    representatives[cluster_key], cluster_spectra.values()))

    logger.info('Export representative evaluation scores to CSV file %s',
                filename_out)
    (pd.DataFrame({'cluster': cluster_keys, **scores})
     .to_csv(filename_out, index=False))


if __name__ == '__main__':
    logging.basicConfig(format='{asctime} [{levelname}/{processName}] '
                               '{module}.{funcName} : {message}',
                        style='{', level=logging.INFO)

    evaluate()

    logging.shutdown()
