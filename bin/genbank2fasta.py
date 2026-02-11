import sys
from Bio import SeqIO


def main(GENBANK: str):
    """ Convert a file in genbank format to FASTA. """
    SeqIO.convert(GENBANK, 'genbank', GENBANK.replace('.gb', '.fa'), 'fasta')
    return

if __name__ == "__main__":
    main(sys.argv[1])
