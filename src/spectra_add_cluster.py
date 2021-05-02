import glob
import logging
from typing import Tuple

import click
import joblib

import ms_io


logger = logging.getLogger('cluster_representative')


@click.command('spectra_add_cluster',
               help='Export spectra to an MGF file containing cluster '
                    'assignments')
@click.option('--spectra', 'filename_spectra',
              help='Input spectrum file(s) (supported file formats: MGF, mzML,'
                   ' mzXML)',
              required=True)
@click.option('--cluster', 'filename_cluster', nargs=2,
              help='Input cluster assignments and cluster type (supported '
                   'clustering formats: MaRaCluster, spectra-cluster, '
                   'MS-Cluster)',
              required=True)
@click.option('--out', 'filename_out',
              help='Output mzML or MGF file containing the updated spectra '
                   'with associated cluster assignments',
              required=True)
def spectra_add_cluster(filename_spectra: str,
                        filename_cluster: Tuple[str, str],
                        filename_out: str):
    filename_spectra = [fn for pattern in filename_spectra
                        for fn in glob.glob(pattern)]
    logger.info('Read spectra from %d peak file(s)', len(filename_spectra))
    spectra = {
        f'{spectrum.filename}:scan:{spectrum.scan}': spectrum
        for spectrum in joblib.Parallel(n_jobs=-1)(
            joblib.delayed(ms_io.read_spectra_list)(filename)
            for filename in filename_spectra)}

    logger.info('Read clusters from cluster file %s', filename_cluster[0])
    clusters = ms_io.read_clusters(
        filename_cluster[0], filename_cluster[1].lower(), list(spectra.keys()))

    for spectrum in spectra.values():
        spectrum_key = f'{spectrum.filename}:scan:{spectrum.scan}'
        if spectrum_key in clusters:
            spectrum.cluster = clusters[spectrum_key]
        else:
            logger.warning('No cluster assignment found for spectrum %s',
                           spectrum.identifier)
    logging.info('Export the spectra including cluster assignments to spectrum'
                 ' file %s', filename_out)
    # Make sure the exported spectra are grouped by their clusters.
    ms_io.write_spectra(filename_out, sorted(
        spectra.values(), key=lambda spectrum: (spectrum.cluster,
                                                spectrum.identifier)))


if __name__ == '__main__':
    logging.basicConfig(format='{asctime} [{levelname}/{processName}] '
                               '{module}.{funcName} : {message}',
                        style='{', level=logging.INFO)

    spectra_add_cluster()

    logging.shutdown()
