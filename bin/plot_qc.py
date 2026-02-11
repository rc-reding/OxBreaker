#!/usr/bin/env python3

import sys
import os
import csv
import numpy as np
import matplotlib.pyplot as plt


# PARAMS
SAMPLE_ID = 0
LENGTH_ID = 2
DEPTH_ID = 4
COV1X_ID = 7
COV10X_ID = 8


def _sort_data(LABELS: list, LENGTH: list, QUAL: list,
               Q20: list, Q30: list) -> list:
    """
        Sort read length, quality, % of reads >Q20, % of reads >Q30 so
        they correspond to samples in 'LABELS'
    """
    NEW_LABELS = sorted(LABELS)  # Avoid sorting in-place
    NEW_LENGTH = list()
    NEW_QUAL = list()
    NEW_Q20 = list()
    NEW_Q30 = list()
    for SAMPLE in NEW_LABELS:
        IDX = LABELS.index(SAMPLE)
        NEW_LENGTH.append(LENGTH[IDX])
        NEW_QUAL.append(QUAL[IDX])
        NEW_Q20.append(Q20[IDX])
        NEW_Q30.append(Q30[IDX])
    return NEW_LABELS, NEW_LENGTH, np.array(NEW_QUAL), \
        NEW_Q20, NEW_Q30

    
def load_data(PATH: str) -> list:
    """
        Load stats from coverage files
    """
    LABELS = list()
    LENGTH = list()
    QUAL = list()
    Q20 = list()
    Q30 = list()

    QC_F_LIST = [f for f in os.listdir(PATH) if f.find('.txt') > -1]
    for FNAME in QC_F_LIST:
        SAMPLE_LABEL = FNAME.split('_')[0].replace('barcode', '')
        SAMPLE_LABEL = SAMPLE_LABEL.replace('SAUR', '')  #  For MRSA samples
        LABELS.append(SAMPLE_LABEL)
        # Screen file
        for line in open(PATH + FNAME):
            if str('Median read length') in line.strip('\n'):
                LENGTH.append(float(line.strip('\n').split(' ')[-1].replace(',', '')))
            elif str('Median read quality') in line.strip('\n'):
                QUAL.append(float(line.strip('\n').split(' ')[-1]))
            elif str('Q20') in line.strip('\n'):
                Q20.append(float(line.split('(')[1].split(')')[0][:-1]))
            elif str('Q30') in line.strip('\n'):
                Q30.append(float(line.split('(')[1].split(')')[0][:-1]))
    return _sort_data(LABELS, LENGTH, QUAL, Q20, Q30)


def plot_qc(LABELS: list, QUAL: str, Q20: list, Q30: list):
    """
    """
    BAR_WIDTH = 0.95
    ACCURACY = 100 * (1 - 10**-(QUAL/10))  # ERR = 10^-(Q/10)
    fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(12, 5))
    ax.bar(LABELS, Q20, width=BAR_WIDTH,
           color='k', label='>99% accuracy')
    ax.bar(LABELS, Q30, width=BAR_WIDTH,
           color='grey', label='>99.9% accuracy')
    ax2 = ax.twinx()
    ax2.plot(LABELS, ACCURACY, linestyle='None', marker='o', color='r')

    # Format axes
    ax.tick_params(axis='both', direction='out')  # right=True, labelright=True)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_xlim([0-BAR_WIDTH/1.9, len(Q20)-BAR_WIDTH/2])
    ax.legend(frameon=False, loc='upper right', bbox_to_anchor=(0.95, 1.1))

    ax2.tick_params(axis='both', direction='out', labelcolor='r')  # right=True, labelright=True)
    ax2.spines['top'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax2.set_xlim([0-BAR_WIDTH/1.9, len(Q20)-BAR_WIDTH/2])
    ax2.set_ylim([90, 100])

    ax.set_ylabel('% of reads', fontsize=24)
    ax2.set_ylabel('Average accuracy', fontsize=24)
    ax.set_xlabel('Sample No.', fontsize=24)
    fig.tight_layout()
    return fig


# Define PATH
PATH = sys.argv[1]

# Load data
LABELS, LENGTHS, QUAL, Q20, Q30 = load_data(PATH)

# Plot
fig = plot_qc(LABELS, QUAL, Q20, Q30)
fig.savefig(sys.argv[2] + 'samples_qc.pdf',
            format='pdf', bbox_inches='tight')
plt.close(fig)

# Load data (Filtered) -- DATA NO LONGER FILTERED
#LABELS, LENGTHS, QUAL, Q20, Q30 = load_data(PATH + str('filtered/'))
#
## Plot (Filtered)
#fig = plot_qc(LABELS, QUAL, Q20, Q30)
#fig.savefig(sys.argv[2] + 'samples_qc_filtered.pdf',
#            format='pdf', bbox_inches='tight')
#plt.close(fig)
