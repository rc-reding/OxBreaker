#!/usr/bin/env python3

import os
import sys
import gzip as gz
import urllib
from Bio import Entrez
import xml.etree.ElementTree as et


USR_EMAIL = str('oxbreaker@ndm.ox.ac.uk')


def _get_assembly_info(ASMBL_ID: str) -> dict:
    """
        Given an entry ID, retrieve  full summary --- which contains URL
        to fasta sequences.
    """
    summary = Entrez.esummary(db='assembly', id=ASMBL_ID, email=USR_EMAIL)
    asmbl_info = Entrez.read(summary)
    return asmbl_info['DocumentSummarySet']['DocumentSummary'][0]


def _get_assembly_url(asmbl_info: dict) -> str:
    """
        Given a dictionary with all the assembly full report,
        parse FTP path and return link for assembly file in FASTA
        format.
    """
    if asmbl_info['FtpPath_RefSeq'] != str(''):
        remote_dir = asmbl_info['FtpPath_RefSeq']
    else:
        remote_dir = asmbl_info['FtpPath_GenBank']
    remote_filename = os.path.basename(remote_dir) + str('_genomic.fna.gz')
    return os.path.join(remote_dir, remote_filename)


def _get_contig_number(metadata_xml: str) -> int:
    """
        Find how many contigs conform the assembly
    """
    # XML in metadata is incomplete
    xml_string = f"<AssemblyData>{metadata_xml}</AssemblyData>"
    metadata = et.fromstring(xml_string)
    for entry in metadata.find('Stats'):
        if entry.get('category') == str('contig_count'):
            break
    return int(entry.text)


def _get_chr_number(report_url: str) -> int:
    """
        Find how many contigs there are for the chromosome
    """
    report_file = urllib.request.urlretrieve(report_url)[0]
    if os.path.exists(report_file):
        chromosomeN = 0
        for entry in open(report_file, 'r').read().split('\n'):
            if len(entry) > 0:  # EOF == len 0
                if entry[0] != str('#') and str('hromosome') in entry:
                    chromosomeN += 1
                    ASMBL_ID = entry.split('\t')[0]
    return chromosomeN, ASMBL_ID


def get_genbank(fasta: str, CHROMOSOME: bool = False) -> str:
    """
        Given a fasta file with a reference sequence,
        extract access code and retrieve full genbank file.
    """
    if CHROMOSOME is False:
        # Get access code from FASTA header
        with gz.open(fasta, 'r') as FASTA_REF:
            ACCESS_CODE = FASTA_REF.readline().decode().split(' ')[0][1:]
    else:
        # If CHROMOSOME == True, 'fasta' is the ACCESS_CODE
        ACCESS_CODE = fasta
    # Retrieve file
    file = Entrez.efetch(db='nucleotide', rettype='gbwithparts',
                         retmode='text', id=ACCESS_CODE,
                         email=USR_EMAIL)
    genbank_file = open('reference.gb', 'w')
    genbank_file.write(file.read())
    genbank_file.close()
    return str('reference.gb')


def get_fasta(species: str) -> list:
    """
     Use anonymised codes in 'sample_file' to populate NCBI
     and find the corresponding FASTQ filenames.
    """
    if type(species) is list:
        # Candidates already sorted by number of reads mapped:
        # Element 0 has the highest.
        success = False
        for candidate in species:
            # Find assembly ID
            search = Entrez.esearch(db='assembly', term=candidate,
                                    retmax='10', email=USR_EMAIL)
            record = Entrez.read(search)
            for ASMBL_ID in record['IdList']:
                # Retrieve metadata for dataset
                info = _get_assembly_info(ASMBL_ID)
                contigN = _get_contig_number(info['Meta'])
                print(candidate, info, contigN)
                if contigN == 1:
                    remote_file = _get_assembly_url(info)
                    success = True
                    break
                elif str('has-plasmid') in info['PropertyList']:
                    report_url = info['FtpPath_Assembly_rpt']
                    chromosomeN, ASMBL_ID = _get_chr_number(report_url)
                    if chromosomeN == 1:
                        success = True
                        remote_file = report_url
                        break
            if success is True:
                print("Candidate found:", candidate)
                break

        if success is False:
            raise Exception("All candidate assemblies have more than 1 contig.")
    else:
        # Find assembly ID
        search = Entrez.esearch(db='assembly', term=candidate,
                                retmax='10', email=USR_EMAIL)
        record_ID = Entrez.read(search)
        ASMBL_ID = record_ID['IdList'][0]

        # Retrieve metadata for dataset
        info = _get_assembly_info(ASMBL_ID)
        contigN = _get_contig_number(info['Meta'])
        if contigN == 1:
            remote_file = _get_assembly_url(info)
        else:
            raise Exception("Candidate assembly has more than 1 contig.")

    # Download remote file
    if contigN == 1:
        urllib.request.urlretrieve(remote_file, "reference.fna.gz")
        if os.path.exists("reference.fna.gz"):
            get_genbank("reference.fna.gz")
            return True, remote_file
    elif chromosomeN == 1:
        get_genbank(ASMBL_ID, CHROMOSOME=True)
        return True, remote_file
    else:
        return False, remote_file


def parse_kraken2_report(REPORT: str) -> str:
    """
        Parse kraken2 report to extract species with
        highest number of reads mapped onto them.

        NB. Report is ordered from group with highest number
            of reads mapped to lowest.
    """
    PROP_ID = 0
    TYPE_ID = 3
    NAME_ID = 5

    # Pin species
    candidate_species = None
    for entry in open(REPORT, 'r'):
        entry = entry.lstrip(' ').rstrip('\n')  # Sanitize output
        entry = entry.split('\t')

        if entry[TYPE_ID] == str('S') and candidate_species is None:
            candidate_species = entry
        elif entry[TYPE_ID] == str('S') and float(entry[PROP_ID]) > float(candidate_species[PROP_ID]):
            candidate_species = entry

    # Pin construct
    species = candidate_species[NAME_ID].strip(' ')
    print("Species:", species)
    candidate_construct = list()
    prev_candidate = None
    for entry in open(REPORT, 'r'):
        entry = entry.lstrip(' ').rstrip('\n')  # Sanitize output
        entry = entry.split('\t')

        if entry[TYPE_ID] == str('S1') and species in entry[NAME_ID] and\
                len(candidate_construct) == 0:
            candidate_construct.append(entry)
        elif entry[TYPE_ID] == str('S1') and species in entry[NAME_ID] and\
                float(entry[PROP_ID]) > 0:
            candidate_construct.append(entry)
        elif entry[TYPE_ID] == str('S2') and species in entry[NAME_ID]:
            if float(entry[PROP_ID]) > 0.5 * float(prev_candidate[PROP_ID]):
                candidate_construct.insert(len(candidate_construct)-1, entry)
            elif float(entry[PROP_ID]) > 0:
                candidate_construct.append(entry)

        prev_candidate = entry

    return candidate_species[NAME_ID].strip(' '),\
        [ c[NAME_ID].strip(' ') for c in candidate_construct]


def main(REPORT: str):
    """
        Find species with most reads mapped onto them by
        kraken2, and retrieve the assembly file.
    """
    species, construct = parse_kraken2_report(REPORT)
    if len(construct) == 0 or species == construct:
        print("Finding assembly using species name:", species)
        success, link = get_fasta(species)
    else:
        print("Finding assembly using construct(s) name:", construct)
        success, link = get_fasta(construct)

    if success is True:
        print("Assembly retrieved successfully from", link)
    else:
        print("Something went wrong, could not retrieve assembly data from", link)
    return


if __name__ == '__main__':
    # sys.argv[1] == kraken2 report file
    main(sys.argv[1])
