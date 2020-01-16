import collections
from typing import Dict, Iterable

import click
import pandas as pd
import spectrum_utils.spectrum as sus

import ms_io
import logging
import metrics as mx

logger = logging.getLogger('specpride')
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])



@click.command('evaluate_representatives',
    short_help='Compare an mgf with clustered spectra with a representatives mgf by metric x'
)
@click.option(
    '--cluster_members_file_name', '-c',
    help='MGF file containing the cluster member spectra',
    required=True
)
@click.option(
    '--representatives_file_name', '-r',
    help='MGF file containing the representative spectra',
    required=True
)
@click.option(
    '--out_file_name', '-o',
    help='TSV table with distances',
    required=True
)
@click.option(
    '--metric_option', '-m',
    help='Distance metric that will be used',
    default="average_cos_dist"
)
def evaluate_representatives(
    cluster_members_file_name,
    representatives_file_name,
    out_file_name,
    metric_option="average_cos_dist"
):
    if metric_option == "average_cos_dist":
        metric = mx.average_cos_dist
    else:
        raise ValueError("Metric not supported")
    cluster_member_spectra = ms_io.read_cluster_spectra(
        cluster_members_file_name
    )
    clusters = collections.defaultdict(list)
    for cluster_member in cluster_member_spectra.values():
        clusters[cluster_member.cluster].append(cluster_member)
    representative_spectra = ms_io.read_cluster_spectra(
        representatives_file_name,
        usi_present=False
    )
    for cluster, cluster_members in clusters.items():
        representative_spectrum = representative_spectra[cluster]
        metric(representative_spectrum, cluster_members)











@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    pass


cli.add_command(evaluate_representatives)


if __name__ == '__main__':
    try:
        logging.basicConfig(format='{asctime} [{levelname}/{processName}] '
                                   '{module}.{funcName} : {message}',
                            style='{', level=logging.DEBUG, force=True)
    except ValueError:
        # force argument not supported on python 3.6
        logging.basicConfig(format='{asctime} [{levelname}/{processName}] '
                                   '{module}.{funcName} : {message}',
                            style='{', level=logging.DEBUG)

    cli()

    logging.shutdown()
