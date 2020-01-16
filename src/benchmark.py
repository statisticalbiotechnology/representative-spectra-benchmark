import collections
# from typing import Dict, Iterable
import click
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
    help='JSON file name that saves the distances',
    required=True
)
@click.option(
    '--metric_option', '-m',
    help='Distance metric that will be used',
    default="average_cos_dist"
)
@click.option(
    '--verbose', '-v',
    help='Print verbose logging messages',
    default=False
)
def evaluate_representatives(
    cluster_members_file_name,
    representatives_file_name,
    out_file_name,
    metric_option="average_cos_dist",
    verbose=False,
):
    if metric_option == "average_cos_dist":
        metric = mx.average_cos_dist
    else:
        raise ValueError("Metric not supported")

    logging.info(f"Reading cluster member spectra from {cluster_members_file_name}")
    cluster_member_spectra = ms_io.read_cluster_spectra(
        cluster_members_file_name
    )
    clusters = collections.defaultdict(list)
    logging.info(f"Grouping cluster member spectra")
    for cluster_member in cluster_member_spectra.values():
        clusters[cluster_member.cluster].append(cluster_member)
    logging.info(f"Reading representatative member spectra from {representatives_file_name}")
    representative_spectra = ms_io.read_cluster_spectra(
        representatives_file_name,
        usi_present=False
    )
    logging.info(f"Calculating distances")
    distances = {}
    i = 0
    for cluster, cluster_members in clusters.items():
        i += 1
        if i > 100:
            break
        if cluster not in representative_spectra:
            if verbose:
                logging.info(f"Cluster: {cluster} has no representative spectrum")
            continue
        representative_spectrum = representative_spectra[cluster]
        distance = metric(representative_spectrum, cluster_members)
        if verbose:
            logging.info(f"Cluster: {cluster}, Distance: {distance}")
        distances[cluster] = distance
    logging.info(f"Saving to JSON file {out_file_name}")
    ms_io.write_distance_dict_to_json(out_file_name, distances)


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
