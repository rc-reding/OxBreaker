#!/usr/bin/env python
import sys
import os
import subprocess
from Bio import SeqIO
from io import StringIO
import numpy as np



def Create_DB(infastapath: str):
    """ Create a blast DB to be used to self-blast assembly. """
    makeblastdb_cmd = 'makeblastdb '
    makeblastdb_args = list(['-dbtype nucl',
                             '-in ' + infastapath,
                             '-out reference_db'])
    os.system(makeblastdb_cmd + ' '.join(makeblastdb_args))
    return 'reference_db'


def Generate_BED(regionspath: str, total_bp: int):
    """
    Produce a BED file based on 'regionspath', so that only
    the sequence outside of repregions are processed downstream
    """
    nonrepeat_fp = open('reference.bed', 'w')
    prev_stop = 1
    prev_repetition = 0
    prev_entry = []
    for region in open(regionspath, 'r'):
        if region[0] != '>':
            rep_start, rep_end = region.strip('\n').split('\t')

            # Skip repetitions if at begining of chromosome, or tandem
            if int(rep_start) != 1 and \
                    prev_stop != int(rep_start) and \
                    prev_stop + 1 != int(rep_start):
                # if repeated region at the begining, skip
                start = int(prev_stop)
                if start == 1:
                    end = int(rep_end) #+ 1
                else:
                    end = int(rep_start) #- 1

                print("Previous repeat:", prev_repetition)
                print("Current Repeat:", rep_start, rep_end)
                print("Non-repeat region:", start, end, '\n')
                assert end <= total_bp
                assert end - start > 0

            # Write data
            if prev_entry == [] or start not in prev_entry:
                nonrepeat_fp.write('{0}\t{1}\t{2}\n'.format(header, str(start), str(end)))
                prev_entry = list([start, end])
            
            # update prev_start/prev_end
            prev_stop = int(rep_end) #+ 1
            prev_repetition = (rep_start, rep_end)
        else:
            header = region[1:].strip('\n')

    # Include end of the chromosome, if relevant
    if prev_stop != total_bp:
        # Write data
        nonrepeat_fp.write('{0}\t{1}\t{2}\n'.format(header, str(prev_stop), str(total_bp)))
    nonrepeat_fp.close()
    return


def Export_Repeats(infastapath: str, regionspath: str):
    """ Export repeated sequences for downstream analysis. Use reference sequence. """
    seq_file = regionspath.replace('.array', '.seq')
    reference = SeqIO.read(infastapath, 'fasta')
    with open(seq_file, 'w') as fOut:
        region_number = 1
        for region in open(regionspath, 'r'):
            if region[0] != '>':
                start, end = region.strip('\n').split('\t')
                fOut.writelines('>Region {0}\n{1}\n'.format(str(region_number), reference.seq[int(start):int(end)]))
                region_number += 1
    return


def Make_Repeat_Mask_Txt(outfastapath: str, prefix: str, word_size=17, window_size=40, gapopen=5, e_thresh=0.0001, perc_identity=90, gapextend=2,
                         min_length=75, threads=10 ):
    """
    Run blastn on contigs in input fasta file against database dbname. Parameters set to NCBI recommended defaults for blastn.
    """
    maskpath = prefix + '_repmask.array'
    regionspath = prefix + '_repregions.array'
    statspath = prefix + '.stats'

    # BLAST cmd line argument
    blast_cmd = list(["blastn",
                      "-query", outfastapath,
                      "-db", prefix,
                      "-evalue", str(e_thresh),
                      "-word_size", str(word_size),
                      "-window_size", str(window_size),
                      "-num_threads", str(threads),
                      "-gapopen", str(gapopen),
                      "-gapextend", str(gapextend),
                      "-dust", "yes",
                      "-perc_identity", str(perc_identity),
                      "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send"
                  ])
    try:
        blast_out, blast_err = subprocess.Popen(blast_cmd, stdout=subprocess.PIPE,
                                                stderr=subprocess.PIPE).communicate()
        assert not blast_err.decode()
    except (AssertionError) as err:
        raise Exception(
    	    'Erro: Blast failed during construction of repeat mask : {0}'.format(err))

    repmask_fp = open(maskpath, 'w')
    repregions_fp = open(regionspath, 'w')
    total_bp = 0
    repetitive_bp = 0
    num_regions = 0

    # each blast_rec is result from one query sequence (contig)
    blast_stream = StringIO(blast_out.decode())
    prev_header = None
    for contig_count, contig in enumerate(SeqIO.parse(outfastapath, 'fasta'), 1):
        if prev_header != contig.name:
            repregions_fp.write('>{0}\n'.format(contig.name))
            prev_header = contig.name
        total_bp += len(contig)
        repmask = np.zeros(len(contig), dtype=bool)
        try:
            fields = next(blast_stream).split()
        except StopIteration:
            fields = None
        while fields and fields[0] == contig.name:
            contig_name, match_name = fields[:2]
            hit_perc_ident = float(fields[2])
            hit_length, q_start, q_end, s_start, s_end = (
                int(x) for x in fields[3:])
            (x1, y1), (x2, y2) = sorted( ((q_start, q_end), tuple(sorted((s_start, s_end)))) )
            if hit_length >= min_length and (contig_name != match_name or not (x2 <= x1 <= y2 and x2 <= y1 <= y2)):
                repmask[q_start - 1:q_end] = True
            try:
                fields = next(blast_stream).split()
            except StopIteration:  # end of blast hits
                fields = None

        # output.bam repmask as 1 and 0, 100 per line
        repmask_fp.write('>{0}\n'.format(contig.name))
        for i in range(0, len(repmask), 100):
            j = min(i + 100, len(repmask))
            repmask_fp.write('{0}\n'.format(''.join(str(i)
                                                    for i in repmask[i:j].astype(int))))
        # identify postitions of repetitive regions (runs of 1s in the
        # repmask array)
        # 0-based numbering
        region_starts = list(np.where(repmask[1:] > repmask[:-1])[0] + 1)
        region_ends = list(np.where(repmask[1:] < repmask[:-1])[0] + 1)
        # special case: full blast hit for this contig against another
        # contig
        if repmask.all():
            region_starts = [0]
            region_ends = [len(repmask)]
        # fix ends, in case regions start from the first position in the
        # sequence or end at the last
        if region_starts and ((not region_ends) or (region_starts[-1] > region_ends[-1])):
            region_ends.append(len(repmask))
        if region_ends and ((not region_starts) or (region_starts[0] > region_ends[0])):
            region_starts = [0] + region_starts

        repregions_fp.writelines('{0}\t{1}\n'.format(
            rs, re) for rs, re in zip(region_starts, region_ends))
        repetitive_bp += repmask.sum()
        num_regions += len(region_starts)

    repmask_fp.close()
    repregions_fp.close()

    #Carlos: Make BED file based on repregions_fp and reference
    Generate_BED(regionspath, total_bp)
    Export_Repeats(outfastapath, regionspath)

    pct_repetitive = '{0:.2f}'.format( (float(repetitive_bp) / total_bp) * 100 )
    statsvalues = '\t'.join((outfastapath, outfastapath, str(contig_count), str(total_bp),\
                             str(repetitive_bp), str(num_regions), pct_repetitive))

    # save result summary
    with open(statspath, 'w') as o:
        o.write('refid\trefcd\tcontigs\tnumbp\trepetitivebp\trepregions\trepetitivepct\n{values}\n'.format(
            values=statsvalues))
    print(
        'Info: Repetitive regions for all of {0}: {1}/{2} bp ({3}%)'.format(outfastapath, repetitive_bp, total_bp,
                                                                            pct_repetitive))
    return


def main(FASTA: str):
    """
    Detect repeated regions in a genome by self=blasting the sequence,
    and produce a bed file containing the non-repeated section of the
    reference
    """
    prefix = Create_DB(FASTA)

    # COMPASS routine defaults to word_size=17, that's too big...?
    # need to see what the repetitions look like
    Make_Repeat_Mask_Txt(FASTA, prefix, word_size=10,
    			 perc_identity=0.95,
    			 window_size=10)
    return


if __name__ == '__main__':
    main(sys.argv[1])
