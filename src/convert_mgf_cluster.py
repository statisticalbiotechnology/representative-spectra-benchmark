import logging
import os

import click
import pyteomics.mgf

import ms_io


logger = logging.getLogger('specpride')


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


def print_help():
    ctx = click.get_current_context()
    click.echo(ctx.get_help())
    ctx.exit()


@click.command('mgf_add_cluster',
               short_help='Add MaRaCluster cluster assignments to MGF')
@click.option('--filename_mgf', '-s', help='MGF file containing the spectra')
@click.option('--filename_cluster', '-c',
              help='File containing MaRaCluster cluster assignments')
@click.option('--filename_out', '-o',
              help='Output MGF file name containing the updated spectra')
@click.option('--px_accession', '-a', help='ProteomeXchange accession of the '
                                           'project (used to compile USIs)')
def mgf_add_cluster(filename_mgf: str, filename_cluster: str,
                    filename_out: str, px_accession: str):
    if (filename_mgf is None or filename_cluster is None or
            filename_out is None or px_accession is None):
        print_help()
        return

    spectra = {}
    for spectrum_dict in pyteomics.mgf.read(filename_mgf):
        spectrum = ms_io._dict_to_spectrum(spectrum_dict)
        spectrum.identifier = ms_io._build_usi(
            px_accession, os.path.basename(os.path.splitext(filename_mgf)[0]),
            spectrum_dict['params']['scans'])
        spectra[spectrum.identifier] = spectrum
    logger.info('Read %d spectra from MGF file %s', len(spectra), filename_mgf)

    clusters = ms_io.read_maracluster_clusters(filename_cluster, px_accession)
    logger.info('Read %d clusters for %d PSMs from cluster file %s',
                clusters.nunique(), len(clusters), filename_cluster)

    logging.info('Export the clustered spectra to MGF file %s', filename_out)
    spectra_dicts = []
    for usi, cluster_i in clusters.items():
        if usi in spectra:
            spectrum = spectra[usi]
            spectrum.identifier = f'cluster-{cluster_i};{spectrum.identifier}'
            spectrum_dict = ms_io._spectrum_to_dict(spectrum)
            spectra_dicts.append(spectrum_dict)
    with open(filename_out, 'w') as f_out:
        pyteomics.mgf.write(spectra_dicts, f_out)


@click.command('convert-mq-marcluster-mzml', short_help='Command to convert MaxQuant Results and MaCluster into MGF')
@click.option('--mq_msms', '-p', help='Peptide information from MaxQuant')
@click.option('--mrcluster_clusters', '-c', help='The information of the clusters from MaRCluster')
@click.option('--mzml_file', '-s', help='The mgf with the corresponding spectra')
@click.option('--output', '-o', help='Output mgf containing the cluster and the spectra information')
@click.option('--px_accession', '-a', help='ProteomeXchange accession of the project')
@click.option('--raw_name', '-r', help='Original name of the RAW file in proteomeXchange')
def convert_mq_mracluster_mzml(mq_msms, mrcluster_clusters, mzml_file, output, px_accession, raw_name):
    raise NotImplementedError   # TODO

    if mq_msms is None or mrcluster_clusters is None or mzml_file is None:
        print_help()

    # Read the input spectra
    exp = MSExperiment()
    MzMLFile().load(mzml_file, exp)
    sl = SpectrumLookup()
    sl.readSpectra(exp, "=(?<SCAN>\\d+)$")

    count = sum(map(lambda spec: spec.getMSLevel() == 2, exp))
    print("Number of Spectra: " + str(count))

    # Read the msms.txt files using, for now the peptides will be a dictionary, where the key is the scan number
    # and the values is the peptide sequence. We need to be aware that we can have cases when one scan can be associated with more
    # than one peptide sequence

    peptides = read_peptides(mq_msms)
    print('Number of Peptides: ' + str(len(peptides)))

    # Read clusters, the clusters will be a map where the key is the scan and the value is the cluster where the scan belongs
    clusters = read_clusters(mrcluster_clusters)
    print("Number of Clusters: " + str(len(clusters)))

    new_exp = MSExperiment(exp)
    new_exp.clear(False)
    for scan in clusters:
        print('scan: ' + str(scan))
        index = sl.findByScanNumber(scan)
        spectra = MSSpectrum(exp[index])
        cluster_accession = clusters[scan]
        if scan in peptides:
            peptide_sequence = peptides[scan]
            spectra.setMetaValue("Peptide sequence", peptide_sequence.encode("UTF-8"))
        spectra.setMetaValue("Cluster accession", cluster_accession.encode("UTF-8"))

        new_exp.addSpectrum(spectra)

    MzMLFile().store(output, new_exp)


@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    pass


cli.add_command(mgf_add_cluster)
# cli.add_command(convert_mq_mracluster_mzml)

if __name__ == '__main__':
    logging.basicConfig(format='{asctime} [{levelname}/{processName}] '
                               '{module}.{funcName} : {message}',
                        style='{', level=logging.DEBUG, force=True)

    cli()

    logging.shutdown()
