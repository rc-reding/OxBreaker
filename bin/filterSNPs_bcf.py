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
    return bool(len(tmp_record[REF_ID]) == len(tmp_record[ALT_ID]) and \
    		str('.') not in tmp_record[ALT_ID])


def is_multiallelic(record: str) -> bool:
    """
      Determine whether a record has MNVs based on the length of the
      sequence at a given position in the reference genome.
    """
    ALT_ID = 4
    META_ID = 7
    DEPTH_ID = 0
    COUNT_ID = -2
    
    tmp_record = record.decode().strip('\n').split('\t')
    read_metadata = tmp_record[META_ID].split(';')
    if str('.') == tmp_record[ALT_ID]:
        if str('0') == read_metadata[COUNT_ID].split(',')[-2] or \
            str('0') == read_metadata[COUNT_ID].split(',')[-1]:
            return False
        else:
            DEPTH = int(read_metadata[DEPTH_ID].split('=')[1])
            ALT_DEPTH = sum([int(c) for c in read_metadata[COUNT_ID].split(',')[-2:]])
            if float(ALT_DEPTH / DEPTH) <= 0.25:
                return False
            else:
                return True
    else:
        return False


def is_gap(record: str) -> bool:
    """
      Determine whether a record is lost based on the total depth
      at a given position in the reference genome.
    """
    META_ID = 7
    DEPTH_ID = 0
    tmp_record = record.decode().strip('\n').split('\t')
    depth_info = tmp_record[META_ID].split(';')
    DEPTH = int(depth_info[DEPTH_ID].split('=')[1])
    return bool(DEPTH == 0)


def is_above_minFreq(record, MIN_FREQ) -> bool:
    """
        Determine whether a record is above the minimum frequency
        set in the pipeline
    """
    META_ID = 7
    DEPTH_ID = 0
    ALT_DEPTH_ID = -2
    tmp_record = record.decode().strip('\n').split('\t')
    depth_info = tmp_record[META_ID].split(';')

    DEPTH = int(depth_info[DEPTH_ID].split('=')[1])
    ALT_DEPTH = sum([int(read_count) for read_count in depth_info[ALT_DEPTH_ID].split(',')[-2:]])
    
    af_freq = float(ALT_DEPTH/DEPTH)
    return bool(af_freq >= float(MIN_FREQ))


def is_above_minDepth(record, MIN_DEPTH) -> bool:
    """
        Determine whether a record is above the minimum read
        support set in the pipeline
    """
    META_ID = 7
    DEPTH_ID = 0
    ALT_DEPTH_ID = -2
    tmp_record = record.decode().strip('\n').split('\t')
    depth_info = tmp_record[META_ID].split(';')
    DEPTH = int(depth_info[DEPTH_ID].split('=')[1])
    if DEPTH < MIN_DEPTH:
        return False
        
    ALT_DEPTH = sum([int(read_count) for read_count in depth_info[ALT_DEPTH_ID].split(',')[-2:]])
    return bool(ALT_DEPTH >= float(MIN_DEPTH))


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
    var_location = int(tmp_record[LOC_ID])
    ALLOWED_COORD = False
    for LOC in NON_REPEATED_COORDINATES:
        if LOC[0] <= var_location <= LOC[1]:
            ALLOWED_COORD = True
    return ALLOWED_COORD


def mask_position(record, gap=False) -> str:
    """
        Mask variant to use for consensus generation,
        rather than calling the reference.
    """
    VARIANT_ID = 4
    tmp_record = record.decode().split('\t')
    if gap is True:
        tmp_record[VARIANT_ID] = str('-')
    else:
        tmp_record[VARIANT_ID] = str('N')
    record = str('\t').join(tmp_record).encode()
    return record


def main(vcf_file: str, bed_file: str, min_freq: str, min_depth: str):
    fOut = vcf_file.replace('.vcf.gz', '_filtered.vcf.gz')
    indel = 0
    depth = 0
    freq = 0
    multial = 0
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
                    # INDEL?
                    if str('INDEL') not in record.decode():
                        if is_snp(record):
                            if not is_above_minDepth(record, min_depth):
                                record = mask_position(record)
                                raw.write(record)
                                depth += 1
                            if not is_above_minFreq(record, min_freq):
                                record = mask_position(record)
                                raw.write(record)
                                freq += 1
                            elif is_above_minFreq(record, min_freq) and \
                                 is_above_minDepth(record, min_depth):
                                raw.write(record)
                        elif is_multiallelic(record):
                            record = mask_position(record)
                            raw.write(record)
                            multial += 1
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
                elif str('INDEL') not in record.decode():
                    # Check indels
                    if is_snp(record) and \
                        is_above_minFreq(record, min_freq) and \
                        is_above_minDepth(record, min_depth) and \
                        check_location(record, non_repet_coords):
                        filtered.write(record)

    print("[" + str(indel+depth+freq) + " variants filtered] '" + str(indel) + "' INDELs, '" + str(depth) +
          "' variants below MIN_DEPTH, '" + str(freq) + "' variants below MIN_FREQ, and '" + str(multial) +
          "' ambiguous multi-allelic sites.\n")
    return


if __name__ == '__main__':
    # argv[1]: reference vcf.gz sequence
    # argv[2]: BED file containing non-repetitive regions
    # argv[3]: min allele frequency
    # argv[4]: min allele frequency
    main(sys.argv[1], sys.argv[2], float(sys.argv[3]), int(sys.argv[4]))
