#!/usr/bin/env python3
from matplotlib import pyplot as plt
import numpy as np
import gzip as gz
import os
import sys

VCF_PATH = str('/vcf/')


def get_snp_depth(VCF_FILE: str) -> list:
    """
        Retrieve the depth of each SNP, and return
        a numpy array
    """
    with gz.open(VCF_FILE, 'rb') as VCF_DATA:
        snp_depth = list()
        for entry in VCF_DATA:
            if entry.decode()[0] != str('#'):
                metadata = entry.decode().strip('\n').split('\t')[-1]
                variant_depth = metadata.split(':')[3].split(',')[1]
                snp_depth.append(int(variant_depth))
    return np.array(snp_depth)


def get_snp_qual(VCF_FILE: str) -> list:
    """
        Retrieve the depth of each SNP, and return
        a numpy array
    """
    with gz.open(VCF_FILE, 'rb') as VCF_DATA:
        snp_depth = list()
        for entry in VCF_DATA:
            if entry.decode()[0] != str('#'):
                metadata = entry.decode().strip('\n').split('\t')[-1]
                variant_depth = metadata.split(':')[1]
                snp_depth.append(int(variant_depth))
    return np.array(snp_depth)


def scatter_plot(PATH: str):
    """
        Produce a scatter plot of Depth vs Quality
    """
    # Sanitize input
    if PATH[-1] != str('/'):
        PATH += str('/')
    # Create figure
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8,7))

    # Iterate through each directory, excluding subsampled data,
    # and load corresponding vcf data
    species = [folder for folder in os.listdir(PATH) if 'DS_Store' not in folder]
    SP_COLOURS = plt.get_cmap('tab20c', len(species)).colors
    for i, sps in enumerate(species):
        pattern = sps[:4]
        for file in [f for f in os.listdir(PATH + sps + VCF_PATH) if pattern in f and 'DS_Store' not in f]:
            if str('x.') not in file:  # exclude subsampled data
                # Define vcf file
                VCF_FILE = PATH + sps + VCF_PATH + file

                sps_initial = sps[0] + str('. ')
                sps_name = sps.split('_')[0][1:]

                # Load snp data
                SNP_DEPTH = get_snp_depth(VCF_FILE)
                SNP_QUAL = get_snp_qual(VCF_FILE)

                ax.scatter(SNP_DEPTH, SNP_QUAL, color=SP_COLOURS[i],
                           alpha=0.85, s=20, label=sps_initial + sps_name)

    # Set labels
    ax.set_ylabel('Variant Quality (Phred-Score)', fontsize=20)
    ax.set_xlabel('Variant Depth (n. reads)', fontsize=20)
    ax.tick_params(axis='both', which='major', labelsize=16)

    # Set legend
    lgn = ax.legend(loc='lower right', frameon=False)

    
    plt.show()
    return


if __name__ == '__main__':
    scatter_plot(sys.argv[1])
