#!/usr/bin/env python3

import gzip
import sys


def is_snp(record: str) -> bool:
    """
      Determine whether a record is a SNP based on the length of the
      sequence at a given position in the reference genome.
    """
    REF_ID = 3
    ALT_ID = 4
    tmp_record = record.decode().split('\t')
    return bool(len(tmp_record[REF_ID]) == len(tmp_record[ALT_ID]))


def is_pass(record):
    """
    Determine if current VCF record is PASS / P, or not, depending upon constraints given
    by variant caller
    """
    PASS_ID = 6
    PASS = False
    tmp_record = record.decode().split('\t')
    if tmp_record[PASS_ID] == str('PASS') or tmp_record[PASS_ID] == str('P'):
        if str('RefCall') not in record.decode():
            PASS = True
    return PASS


def is_pileup(record) -> bool:
    """
    Determine if current VCF record comes from the alignment (F) or pileup (P).
    !!!!!_ONLY WORKS WITH CLAIR3_!!!!
    """
    PILEUP_ID = 7
    tmp_record = record.decode().split('\t')
    return bool(tmp_record[PILEUP_ID] == str('P'))


def is_above_minFreq(record, MIN_FREQ) -> bool:
    """
        Determine whether a record is above the minimum frequency
        set in the pipeline
    """
    META_ID = -1
    FREQ_ID = -1
    tmp_record = record.decode().strip('\n').split('\t')
    freq_info = tmp_record[META_ID].split(':')
    af_freq = float(freq_info[FREQ_ID])
    return bool(af_freq >= float(MIN_FREQ))


def load_non_rep_locations(BED_FILE) -> list:
    """
        Load coordinates from BED file containing the coordinates
        analysed.
    """
    NON_REPEATED_COORDINATES = list()
    for entry in open(BED_FILE, 'r'):
        _, init, end = entry.strip('\n').split('\t')
        NON_REPEATED_COORDINATES.append([int(init), int(end)])
    return NON_REPEATED_COORDINATES


def check_location(record, NON_REPEATED_COORDINATES) -> bool:
    """
        Determine whether a record is located in a non-repetitive
        region. If not, return False.
    """
    LOC_ID = 1
    tmp_record = record.decode().strip('\n').split('\t')
    var_location = int(tmp_record[META_ID])
    ALLOWED_COORD = False
    for LOC in NON_REPEATED_COORDINATES:
        if LOC[0] <= var_location <= LOC[1]:
            ALLOWED_COORD = True
    return ALLOWED_COORD


def mask_position(record) -> str:
    """
        Mask variant to use for consensus generation,
        rather than calling the reference.
    """
    VARIANT_ID = 4
    tmp_record = record.decode().split('\t')
    tmp_record[VARIANT_ID] = str('N')
    record = str('\t').join(tmp_record).encode()
    return record


def main(vcf_file: str, bed_file: str, min_freq: str, min_depth: str):
    fOut = vcf_file.replace('.vcf.gz', '_filtered.vcf.gz')
    indel = 0
    align = 0
    freq = 0
    # Load coordinates analysed by the pipeline
    non_repet_coords = load_non_rep_locations(bed_file)

    # Report of _all_ variants found, masked where appropriate
    # to use for consensus sequence
    with gzip.open(fOut, 'w') as raw:
        with gzip.open(vcf_file, 'r') as variants:
            for record in variants:
                if record.decode()[0] == str('#'):
                    # Write header as-is
                    raw.write(record)
                elif check_location(record, non_repet_coords):
                    # Check indels
                    if is_snp(record) and is_pileup(record) and \
                            is_above_minFreq(record, min_freq) and is_pass(record):
                        raw.write(record)
                    elif is_snp(record) and is_pileup(record) and \
                            not is_above_minFreq(record, min_freq) and is_pass(record):
                        freq += 1
                        record = mask_position(record)
                        raw.write(record)
                    elif is_snp(record) and not is_pileup(record):
                        align += 1
                    else:
                        indel += 1

    # Report of only sequences compliant with constraints
    # to use for report and genomic context generation
    with gzip.open(fOut.replace('.vcf.gz', '_report.vcf.gz'), 'w') as filtered:
        with gzip.open(vcf_file, 'r') as variants:
            for record in variants:
                if record.decode()[0] == str('#'):
                    # Write header as-is
                    filtered.write(record)
                else:
                    # Check indels
                    if is_snp(record) and is_pileup(record) and \
                            is_above_minFreq(record, min_freq) and is_pass(record):
                        filtered.write(record)

    print("[" + str(indel+align+freq) + " variants filtered] '" + str(indel) + "' INDELs, '" + str(align) +
          "' variants found by alignment only, and '" + str(freq) + "' variants below MIN_FREQ.\n")
    return


if __name__ == '__main__':
    # argv[1]: reference vcf.gz sequence
    # argv[2]: BED file containing non-repetitive regions
    # argv[3]: min allele frequency
    # argv[4]: min allele frequency
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
