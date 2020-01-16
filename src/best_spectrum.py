import collections
from typing import Dict, Iterable

import click
import pandas as pd
import spectrum_utils.spectrum as sus

import ms_io
import logging

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


def get_best_representative(spectra: Dict[str, sus.MsmsSpectrum],
                            scores: pd.Series) -> sus.MsmsSpectrum:
    """
    Get the best representative spectrum with the highest score among the given
    spectra.

    If none of the given spectra have a score a ValueError is raised.

    If two spectra have the same maximum score, the spectrum with the
    alphanumeric first file name with the lowest scan number is selected.

    Parameters
    ----------
    spectra : Dict[str, sus.MsmsSpectrum]
        A dictionary with as keys the USIs and as values the corresponding
        spectra in a cluster.
    scores : pd.Series
        A `Series` with as index the USI and as values the identification
        scores.

    Returns
    -------
    sus.MsmsSpectrum
        The representative spectrum with the highest score.

    Raises
    ------
    ValueError
        If none of the given spectra have a score.
    """
    scores = scores[scores.index.isin(spectra.keys())]
    if len(scores) == 0:
        raise ValueError('No scores found for the given scan numbers')
    return spectra[scores.idxmax()]


def split_into_clusters(spectra: Dict[str, sus.MsmsSpectrum]) \
        -> Iterable[Dict[str, sus.MsmsSpectrum]]:
    """
    Yield collections of spectra grouped by cluster membership.

    Parameters
    ----------
    spectra : Dict[str, sus.MsmsSpectrum]
        The spectra to be grouped by cluster membership. The `MsmsSpectrum`
        objects should have a member variable `cluster` to designate which
        cluster they belong to.

    Returns
    -------
    Iterable[Dict[str, sus.MsmsSpectrum]]
        A dictionary with as keys the USI and as values the spectra for each
        individual cluster.
    """
    clusters = collections.defaultdict(list)
    for spectrum in spectra.values():
        clusters[spectrum.cluster].append(spectrum.identifier)
    for cluster_members in clusters.values():
        yield {usi: spectra[usi] for usi in cluster_members}


@click.command('best_spectrum',
               short_help='Select cluster representatives with the highest '
                          'search engine score')
@click.option('--filename_mgf_in', '-s',
              help='MGF file containing the spectra')
@click.option('--filename_mgf_out', '-o',
              help='Output MGF file name containing the cluster '
                   'representatives')
@click.option('--filename_msms', '-m',
              help='File containing MaxQuant identifications (msms.txt)')
@click.option('--px_accession', '-a', help='ProteomeXchange accession of the '
                                           'project (used to compile USIs)')
def best_spectrum(filename_mgf_in: str, filename_mgf_out: str,
                  filename_msms: str, px_accession: str) -> None:
    """
    Represent clusters by their highest scoring member.

    Parameters
    ----------
    filename_mgf_in : str
        The file name of the MGF input file that contains the original spectra.
    filename_mgf_out : str
        The file name of the MGF output file that contains the cluster
        representatives.
    filename_msms : str
        The file name of the MaxQuant identifications file containing PSM
        scores.
    px_accession : str
        ProteomeXchange accession of the project (used to compile USIs).
    """
    scores = ms_io.read_maxquant_psms(filename_msms, px_accession)['Score']
    spectra = ms_io.read_cluster_spectra(filename_mgf_in)
    representatives = []
    for cluster in split_into_clusters(spectra):
        try:
            representative = get_best_representative(cluster, scores)
            representative.identifier = (f'{representative.cluster};'
                                         f'{representative.identifier}')
            representatives.append(ms_io._spectrum_to_dict(representative))
        except ValueError:
            logging.warning(f"Failed to find best spectrum for cluster: {cluster}")
    ms_io.write_mgf(filename_mgf_out, representatives)


@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    pass


cli.add_command(best_spectrum)
# cli.add_command(convert_mq_mracluster_mzml)

if __name__ == '__main__':
        cli()