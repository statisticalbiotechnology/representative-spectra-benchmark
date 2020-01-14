import click
from pyteomics import mgf

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


def print_help():
    ctx = click.get_current_context()
    click.echo(ctx.get_help())
    ctx.exit()


def buid_usi_accession(cluster_id, peptide_sequence, scan, px_accession, raw_name):
    id = cluster_id + ";" + 'mzspec' + ":" +px_accession + ":" + raw_name + ":" + ":scan:" +str(scan)
    if peptide_sequence is not None:
        id = id + ":" + peptide_sequence
    return id


@click.command('convert-mq-marcluster', short_help='Command to convert MaxQuant Results and MaCluster into MGF')
@click.option('--mq_msms','-p', help='Peptide information from MaxQuant')
@click.option('--mrcluster_clusters','-c', help='The information of the clusters from MaRCluster')
@click.option('--mgf_file', '-s',help='The mgf with the corresponding spectra')
@click.option('--output', '-o', help= 'Output mgf containing the cluster and the spectra information')
@click.option('--px_accession', '-a', help='ProteomeXchange accession of the project')
@click.option('--raw_name', '-r', help='Original name of the RAW file in proteomeXchange')
def convert_mq_mracluster(mq_msms, mrcluster_clusters, mgf_file, output, px_accession, raw_name):

    if(mq_msms is None or mrcluster_clusters is None or mgf_file is None):
        print_help()

    ## Read the input spectra
    input_spectra = mgf.read(mgf_file)
    print('Number of Spectra: ' + str(len(input_spectra)))
    ## Read the msms.txt files using, for now the peptides will be a dictionary, where the key is the scan number
    ## and the values is the peptide sequence. We need to be aware that we can have cases when one scan can be associated with more
    ## than one peptide sequence

    peptides = {}
    with open(mq_msms) as mq_peptides:
        next(mq_peptides)  # skip header
        for line in mq_peptides:
            words = line.split('\t')
            rscan = int(words[1])
            rpept = words[7][1:-1]
            peptides[rscan] = rpept

    print('Number of Peptides: ' + str(len(peptides)))

    ## Read clusters, the clusters will be a map where the key is the scan and the value is the cluster where the scan belongs
    clusters = {}
    cluster_prefix = 'cluster-'
    cluster_index  = 1
    with open(mrcluster_clusters) as cluster_def:
        for line in cluster_def:
            if not line.strip():
                cluster_index = cluster_index + 1
            else:
                words = line.split('\t')
                clusters[int(words[1])] = cluster_prefix + str(cluster_index)
    print("Number of Clusters: " + str(cluster_index))

    for scan in clusters:
        print('scan: ' + str(scan))
        for spectra in input_spectra:
            print(spectra)
            if str(scan) in spectra['params']['title']:
                cluster_accession = clusters[scan]

                if scan not in peptides:
                    peptide_sequence = None
                else:
                    peptide_sequence = peptides[scan]
                spectra['params']['title'] = buid_usi_accession(cluster_accession, peptide_sequence, scan, px_accession, raw_name)
                mgf.write([spectra], output)

@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    """This is the main tool that give access to all commands and options provided by the pypgatk"""


cli.add_command(convert_mq_mracluster)

if __name__ == "__main__":
    cli()


