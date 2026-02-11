#!/usr/bin/env python3
from Bio import SeqIO
import sys


def import_genbank(PATH: str) -> list:
    """ Load genbank record. """
    records = SeqIO.read(PATH, 'genbank')
    return records


def calculate_CDS(RECORDS) -> int:
    """ Calculate length, in nt, of all coding sequences combined. """
    LEN_GENOME = len(RECORDS)
    LEN_CDS = 0

    for entry in RECORDS.features:
        if entry.type == str('CDS'):
            if entry.location.parts == 1:
                LEN_CDS += entry.location.end - entry.location.start
            else:
                for part in entry.location.parts:
                    LEN_CDS += part.end - part.start
    return LEN_GENOME, LEN_CDS


def main(PATH: str):
    """

    """
    records = import_genbank(PATH)
    genome_nt, cds_nt = calculate_CDS(records)
    with open('reference.stats', 'w') as fOut:
        fOut.write('Genome size (nt)\tCoding Genome (nt)\tFraction Coding\tFraction Non-Coding\n')
        fOut.write('{0}\t{1}\t{2}\t{3}\n'.format(str(genome_nt), str(cds_nt),
                                               str(cds_nt / genome_nt),
                                               str(1 - (cds_nt / genome_nt))))
    return


if __name__ == '__main__':
    main(sys.argv[1])
