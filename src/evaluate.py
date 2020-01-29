import logging

import click
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
              help='Output JSON file containing representative scores')
@click.option('--measure', 'measure', type=click.Choice(['avg_dot']),
              help='Measure used to evaluate cluster representatives (options:'
                   ' "avg_dot")')
def evaluate(filename_spectra: str, filename_representatives: str,
             filename_out: str, metric: str = 'avg_dot') -> None:
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
    """
    if metric == 'avg_dot':
        metric_func = metrics.average_cos_dist
    else:
        raise ValueError('Unknown metric specified: %s', metric)

    logger.info('Read original spectra from spectrum file %s',
                filename_spectra)
    spectra = {f'{spectrum.filename}:scan:{spectrum.scan}': spectrum
               for spectrum in ms_io.read_spectra(filename_spectra)}

    logger.info('Read cluster representatives from spectrum file %s',
                filename_representatives)
    representatives = {
        spectrum.cluster: spectrum
        for spectrum in ms_io.read_spectra(filename_representatives)}

    logger.info('Evaluate cluster representatives')
    representative_scores = {}
    for cluster_key, cluster_spectra in tqdm.tqdm(
            representative.get_cluster_spectra(spectra),
            desc='Clusters processed', unit='clusters'):
        if cluster_key in representatives:
            representative_scores[cluster_key] = metric_func(
                representatives[cluster_key], cluster_spectra.values())

    logger.info('Export representative evaluation scores to JSON file %s',
                filename_out)
    ms_io.write_json(filename_out, representative_scores)


if __name__ == '__main__':
    logging.basicConfig(format='{asctime} [{levelname}/{processName}] '
                               '{module}.{funcName} : {message}',
                        style='{', level=logging.INFO)

    evaluate()

    logging.shutdown()
