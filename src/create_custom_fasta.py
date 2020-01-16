from pyteomics import parser, fasta
import argparse

def smart_reverse(prot_seq):
    '''
    Produces a reverse decoy sequence keeping all K and R aminoacids in the same positions.
    '''
    return ''.join(fasta.decoy_sequence(pep, mode='reverse', keep_cterm=True) for pep in parser._cleave(prot_seq, '[RK]', 0))

def make_reverse_fasta(input_file, output_file):
    '''
    Takes as input fasta file, drops all _REVERSED proteins and creates a new _REVERSED decoy proteins.
    '''
    prots = []
    for prot_desc, prot_seq in fasta.read(input_file):
        if not prot_desc.endswith('_REVERSED'):
            prots.append((prot_desc, prot_seq))
            prots.append((prot_desc+'_REVERSED', smart_reverse(prot_seq)))
    fasta.write(prots, output_file, file_mode='w')

def main():
    pars = argparse.ArgumentParser()
    pars.add_argument('input', help='input fasta')
    pars.add_argument('output', nargs='?', help='Output fasta file (default is stdout).')
    args = pars.parse_args()
    make_reverse_fasta(args.input, args.output)

if __name__ == '__main__':
    main()