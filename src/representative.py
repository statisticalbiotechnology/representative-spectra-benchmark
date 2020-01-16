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
    clusters = collections.defaultdict(list)
    for spectrum_key, spectrum in spectra.items():
        clusters[spectrum.cluster].append(spectrum_key)
    for cluster_key, cluster_members in clusters.items():
        yield cluster_key, {spectrum_key: spectra[spectrum_key]
                            for spectrum_key in cluster_members}


@click.command('representative',
               help='Export representative spectra for each cluster to an MGF '
                    'file')
@click.option('--in', '-i', 'filename_in',
              help='Input MGF file containing cluster assignments',
              required=True)
@click.option('--out', '-o', 'filename_out',
              help='Output MGF file containing representative spectra for each'
                   ' cluster',
              required=True)
@click.option('--representative_method', '-r', 'representative_method',
              help='Method used to select the representative spectrum for each'
                   ' cluster (options: "most_similar", "best_spectrum")',
              required=True,
              type=click.Choice(['most_similar', 'best_spectrum']))
@click.option('--distance', '-d', 'distance',
              help='Distance metric to compare spectra to each other (options:'
                   ' "dot")',
              type=click.Choice(['dot']))
@click.option('--psm', '-p', 'filename_psm',
              help='Input PSM file (optional; required for specific '
                   'representative selectors')
def representative(filename_in: str, filename_out: str,
                   representative_method: str, distance: str = None,
                   filename_psm: str = None) \
        -> None:
    logger.info('Read spectra from spectrum file %s', filename_in)
    spectra = {f'{spectrum.filename}:scan:{spectrum.scan}': spectrum
               for spectrum in ms_io.read_spectra(filename_in)}

    # TODO: Other selectors.
    if representative_method == 'most_similar':
        rs = selector.MostSimilarRepresentativeSelector(distance)
        logger.info('Select cluster representatives using the most similar '
                    'spectrum')
    elif representative_method == 'best_spectrum':
        rs = selector.BestSpectrumRepresentativeSelector(filename_psm)
        logger.info('Select cluster representatives using the best spectrum')
    else:
        raise ValueError('Unknown method to select the representative spectra')

    cluster_representatives = []
    for cluster_key, cluster_spectra in tqdm.tqdm(
            get_cluster_spectra(spectra), desc='Clusters processed',
            unit='clusters'):
        try:
            cluster_representative = rs.select_representative(cluster_spectra)
            cluster_representative.identifier = (f'{rs.get_description()}:'
                                                 f'cluster={cluster_key}')
            cluster_representatives.append(cluster_representative)
        except ValueError:
            # logger.warning('No representative found for cluster %s',
            #                cluster_key)
            pass

    logger.info('Export %d cluster representatives to MGF file %s',
                len(cluster_representatives), filename_out)
    ms_io.write_spectra(filename_out, cluster_representatives)


if __name__ == '__main__':
    logging.basicConfig(format='{asctime} [{levelname}/{processName}] '
                               '{module}.{funcName} : {message}',
                        style='{', level=logging.DEBUG)

    representative()

    logging.shutdown()
