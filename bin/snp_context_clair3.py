#!/usr/bin/env python3

import gzip
import sys
from Bio import SeqIO


LOC_ID = 1
REF_ID = 3
ALT_ID = 4
META_ID = 9


def load_genbank(genbank_file: str) -> list:
    """ Load genbank metadata, return CDS coordinates. """
    records = SeqIO.read(genbank_file, 'genbank')
    coordinates = list()
    for FEAT in records.features:
        if FEAT.type == str('CDS'):
            coordinates.append([FEAT.location.start, FEAT.location.end])
    return coordinates


def is_cds(location: int, cds_coordinates: list) -> bool:
    """ Check if location is within the cds coordinates. """
    IS_CDS = [True for c in cds_coordinates if c[0] <= location <= c[1]]
    if len(IS_CDS) == 0:
        return False
    else:
        return True


def main(vcf_file: str, reference_file: str, genbank_file: str):
    """
    Extract ±4nt either side of the SNP to know its context
    and analyse downstream possible repetitions
    """
    fOut = vcf_file.replace('.vcf.gz', '_variants_context.tsv')
    reference = SeqIO.read(reference_file, 'fasta')
    cds_coordinates = load_genbank(genbank_file)

    with open(fOut, 'w') as ctx_fOut:
        # Write header
        ctx_fOut.write('Range\tSubstitution\tContext\tLocation\tFrequency\n')

        with gzip.open(vcf_file, 'r') as variants:
            ctx_number = 1
            for record in variants:
                if record.decode()[0] != str('#'):
                    SNP_LOCATION = int(record.decode().split('\t')[LOC_ID])
                    IS_CDS = is_cds(SNP_LOCATION, cds_coordinates)

                    LOCATION_LABEL = str('CDS') if IS_CDS is True else str('Intergenic')

                    REF_NT = record.decode().split('\t')[REF_ID]
                    ALT_NT = record.decode().split('\t')[ALT_ID]
                    ctx_seq = reference.seq[SNP_LOCATION-5:SNP_LOCATION+4]

                    METADATA = record.decode().strip('\n').split('\t')[META_ID]
                    FREQ = float(METADATA.split(':')[-1]) * 100

                    # Write context
                    ctx_fOut.write('{0}-{1}\t{2}->{3}\t{4}\t{5}\t{6}\n'.format(str(SNP_LOCATION - 4),
                                                                                str(SNP_LOCATION + 5),
                                                                                REF_NT, 
                                                                                ALT_NT, 
                                                                                ctx_seq,
                                                                                LOCATION_LABEL,
                                                                                str(FREQ)))
                    ctx_number += 1
    return


if __name__ == '__main__':
    # argv[1]: reference vcf.gz sequence
    # argv[2]: min allele frequency
    main(sys.argv[1], sys.argv[2], sys.argv[3])
