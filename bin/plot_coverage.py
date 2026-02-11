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
COV1X_ID = 8
COV10X_ID = 9


def _sort_data(LABELS: list, LENGTH: list, DEPTH: list,
               COV1X: list, COV10X: list) -> list:
    """
        Sort length, depth, and breadth of coverage so
        they correspond to samples in 'LABELS'
    """
    NEW_LABELS = sorted(LABELS)  # Avoid sorting in-place
    NEW_LENGTH = list()
    NEW_DEPTH = list()
    NEW_COV1X = list()
    NEW_COV10X = list()
    for SAMPLE in NEW_LABELS:
        IDX = LABELS.index(SAMPLE)
        NEW_LENGTH.append(LENGTH[IDX])
        NEW_DEPTH.append(DEPTH[IDX])
        NEW_COV1X.append(COV1X[IDX])
        NEW_COV10X.append(COV10X[IDX])
    return NEW_LABELS, NEW_LENGTH, NEW_DEPTH, \
        NEW_COV1X, NEW_COV10X


def load_data(PATH: str) -> list:
    """
        Load stats from coverage files
    """
    LABELS = list()
    LENGTH = list()
    DEPTH = list()
    COV1X = list()
    COV10X = list()

    DEPTH_FILE_LIST = [f for f in os.listdir(PATH) if f.find('.csv') > -1]
    for FNAME in DEPTH_FILE_LIST:
        print(PATH + FNAME)
        raw_data = csv.reader(open(PATH + FNAME, 'r'))
        next(raw_data)  # Skip header
        raw_data = next(raw_data)  # Hack since 'raw_data' is an Iterator
        SAMPLE_LABEL = raw_data[SAMPLE_ID].replace('barcode', '')
        SAMPLE_LABEL = SAMPLE_LABEL.replace('SAUR', '')  # MRSA
        LABELS.append(SAMPLE_LABEL)
        LENGTH.append(int(raw_data[LENGTH_ID]))
        DEPTH.append(float(raw_data[DEPTH_ID]))
        COV1X.append(float(raw_data[COV1X_ID]) if len(raw_data[COV1X_ID]) > 0 else 0.0)
        COV10X.append(float(raw_data[COV10X_ID]) if len(raw_data[COV10X_ID]) > 0 else 0.0)
    return _sort_data(LABELS, LENGTH, DEPTH, COV1X, COV10X)


def plot_breadth_coverage(LABELS: list, COV1X: list, COV10X: list):
    """
    """
    BAR_WIDTH = 0.95
    fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(12, 5))
    ax.bar(LABELS, np.multiply(COV1X, 100), width=BAR_WIDTH,
           color='k', label='At least 1X')
    ax.bar(LABELS, np.multiply(COV10X, 100), width=BAR_WIDTH,
           color='grey', label='At least 10X')
    ax.plot([0-BAR_WIDTH/2, len(COV1X)-BAR_WIDTH/2], [100, 100], linestyle='--',
            linewidth=0.5, color='k')

    # Format axes
    ax.tick_params(axis='both', direction='out', right=True, labelright=True)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_xlim([0-BAR_WIDTH/1.9, len(COV1X)-BAR_WIDTH/2])
    ax.legend(frameon=False, loc='upper right', bbox_to_anchor=(0.95, 1.1))

    plt.ylabel('Reference covered ($\%$)', fontsize=24)
    plt.xlabel('Sample No.', fontsize=24)
    plt.tight_layout()
    return fig


def plot_depth(LABELS: list, DEPTHS: list):
    """
    """
    BAR_WIDTH = 0.95
    fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(12, 5))
    ax.bar(LABELS, DEPTHS, width=BAR_WIDTH,
           color='k')

    # Format axes
    ax.tick_params(axis='both', direction='out', right=True, labelright=True)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_xlim([0-BAR_WIDTH/1.9, len(COV1X)-BAR_WIDTH/2])
    ax.legend(frameon=False, loc='upper right', bbox_to_anchor=(0.95, 1.1))

    plt.ylabel('Average Coverage Depth', fontsize=24)
    plt.xlabel('Sample No.', fontsize=24)
    plt.tight_layout()
    return fig

# Define PATH
PATH = sys.argv[1]

# Load data
LABELS, _, DEPTHS, COV1X, COV10X = load_data(PATH)

# Plot
fig = plot_breadth_coverage(LABELS, COV1X, COV10X)
fig.savefig(sys.argv[2] + 'samples_breadth_coverage.pdf',
            format='pdf', bbox_inches='tight')
plt.close(fig)

fig = plot_depth(LABELS, DEPTHS)
fig.savefig(sys.argv[2] + 'samples_coverage_depth.pdf',
            format='pdf', bbox_inches='tight')
plt.close(fig)

