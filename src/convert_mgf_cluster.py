import re

import click
from pyteomics import mgf, mzml
from pyopenms import *

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


def print_help():
    ctx = click.get_current_context()
    click.echo(ctx.get_help())
    ctx.exit()


def buid_usi_accession(cluster_id, peptide_sequence, scan, px_accession, raw_name, charge):
    usi = cluster_id + ";" + 'mzspec' + ":" + px_accession + ":" + raw_name + ":" + "scan:" + str(scan)
    if peptide_sequence is not None:
        usi = usi + ":" + peptide_sequence + "/" + str(charge)
    return usi


def read_peptides(mq_msms):
    peptides = {}
    with open(mq_msms) as mq_peptides:
        next(mq_peptides)  # skip header
        for line in mq_peptides:
            words = line.split('\t')
            rscan = int(words[1])
            rpept = words[7][1:-1]
            peptides[rscan] = rpept
    return peptides


def read_clusters(mrcluster_clusters):
    clusters = {}
    cluster_prefix = 'cluster-'
    cluster_index = 1
    with open(mrcluster_clusters) as cluster_def:
        for line in cluster_def:
            if not line.strip():
                cluster_index = cluster_index + 1
            else:
                words = line.split('\t')
                clusters[int(words[1])] = cluster_prefix + str(cluster_index)
    return clusters


@click.command('convert-mq-marcluster', short_help='Command to convert MaxQuant Results and MaCluster into MGF')
@click.option('--mq_msms', '-p', help='Peptide information from MaxQuant')
@click.option('--mrcluster_clusters', '-c', help='The information of the clusters from MaRCluster')
@click.option('--mgf_file', '-s', help='The mgf with the corresponding spectra')
@click.option('--output', '-o', help='Output mgf containing the cluster and the spectra information')
@click.option('--px_accession', '-a', help='ProteomeXchange accession of the project')
@click.option('--raw_name', '-r', help='Original name of the RAW file in proteomeXchange')
def convert_mq_mracluster_mgf(mq_msms, mrcluster_clusters, mgf_file, output, px_accession, raw_name):
    if mq_msms is None or mrcluster_clusters is None or mgf_file is None:
        print_help()

    # Read the input spectra
    input_spectra = mgf.read(mgf_file)
    spectra_list = list(input_spectra)
    spectra_dic = {}
    for spectrum in spectra_list:
        scan = spectrum['params']['title'].split("=")[-1]
        spectra_dic[int(scan)] = spectrum
    print('Number of Spectra: ' + str(len(spectra_list)))

    # Read the msms.txt files using, for now the peptides will be a dictionary, where the key is the scan number
    # and the values is the peptide sequence. We need to be aware that we can have cases when one scan can be associated with more
    # than one peptide sequence

    peptides = read_peptides(mq_msms)
    print('Number of Peptides: ' + str(len(peptides)))

    # Read clusters, the clusters will be a map where the key is the scan and the value is the cluster where the scan belongs
    clusters = read_clusters(mrcluster_clusters)
    print("Number of Clusters: " + str(len(clusters)))

    for scan in clusters:
        print('scan: ' + str(scan))
        if scan in spectra_dic:
            cluster_accession = clusters[scan]
            spectrum = spectra_dic[scan]
            if scan not in peptides:
                peptide_sequence = None
            else:
                peptide_sequence = peptides[scan]
            charge = int(spectrum['params']['charge'][0])
            spectrum['params']['title'] = buid_usi_accession(cluster_accession, peptide_sequence, scan, px_accession,
                                                                raw_name, charge)
            mgf.write([spectrum], output)


@click.command('convert-mq-marcluster-mzml', short_help='Command to convert MaxQuant Results and MaCluster into MGF')
@click.option('--mq_msms', '-p', help='Peptide information from MaxQuant')
@click.option('--mrcluster_clusters', '-c', help='The information of the clusters from MaRCluster')
@click.option('--mzml_file', '-s', help='The mgf with the corresponding spectra')
@click.option('--output', '-o', help='Output mgf containing the cluster and the spectra information')
@click.option('--px_accession', '-a', help='ProteomeXchange accession of the project')
@click.option('--raw_name', '-r', help='Original name of the RAW file in proteomeXchange')
def convert_mq_mracluster_mzml(mq_msms, mrcluster_clusters, mzml_file, output, px_accession, raw_name):
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
    """This is the main tool that give access to all commands and options provided by the pypgatk"""


cli.add_command(convert_mq_mracluster_mgf)
cli.add_command(convert_mq_mracluster_mzml)

if __name__ == "__main__":
    cli()

# <cvParam cvRef="MS" accession="MS:1000511" value="1" name="ms level" />
# https://www.ebi.ac.uk/ols/ontologies/ms/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FMS_1000889

# <userParam name="[Thermo Trailer Extra]Monoisotopic M/Z:" type="xsd:float" value="665.8334" />
# <userParam name="Cluster accession" type="xsd:string" value="cluster-1" />
