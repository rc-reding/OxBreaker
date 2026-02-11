#!/usr/bin/env python3
import sys
import os
import gzip as gz
import numpy as np
from matplotlib import pyplot as plt


def get_SNP_data(VCF_FILE: str) -> list:
    """
        Retrieve the depth of each SNP, and return
        a numpy array
    """
    with gz.open(VCF_FILE, 'rb') as VCF_DATA:
        snp_depth = list()
        for entry in VCF_DATA:
            if entry.decode()[0] != str('#'):
                metadata = entry.decode().strip('\n').split('\t')[-1]
                print("Metadata", metadata, entry.decode())
                variant_depth = metadata.split(':')[3].split(',')[1]
                snp_depth.append(int(variant_depth))
    return np.array(snp_depth)


def get_asmbl_depth(ASMBL_METADATA: str) -> float:
    """
        Get mean depth of the assembly and standard deviation
    """
    MEAN_ID = 4
    STD_ID = 5
    with open(ASMBL_METADATA, 'rb') as METADATA:
        METADATA.readline()  # Skip header
        data = METADATA.readline().decode().split(',')
    return float(data[MEAN_ID]), float(data[STD_ID])


def plot_histogram(VCF_FILE: str, SNP_DATA: list, DEPTH: float,
                   DEPTH_STD: float):
    """
        Plot histogram of SNP depth distribution
    """
    # Get relative frequencies of each depth
    hist = np.histogram(SNP_DATA, density=True)
    snp_freq = hist[0] * np.diff(hist[1])  # sum() == 1.0

    # Create figure
    fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(7, 6.5))
    #ax.bar(hist[1][:-1], height=snp_freq, width=1.5,
    #       edgecolor='black', linewidth=1)

    ax.hist(SNP_DATA, width=1.5,
                 edgecolor='black', linewidth=1)
    # Add assembly mean±std
    y_lims = ax.get_ylim()
    ax.plot([DEPTH, DEPTH], y_lims,
            color='red', linewidth=3)
    std_x1_patch = plt.Rectangle(xy=[DEPTH-DEPTH_STD, y_lims[0]],
                              width=2*DEPTH_STD,
                              height=y_lims[1],
                              color='salmon', alpha=0.5,
                              zorder=-5,
                              label='Mean ± 1 x StDev')
    std_x2_patch = plt.Rectangle(xy=[DEPTH-(DEPTH_STD*2), y_lims[0]],
                              width=2*DEPTH_STD*2,
                              height=y_lims[1],
                              color='pink', alpha=0.5,
                              zorder=-10,
                              label='Mean ± 2 x StDev')
    ax.add_patch(std_x1_patch)
    ax.add_patch(std_x2_patch)


    # Restore former y_lim
    ax.set_ylim(y_lims)

    # Legend
    lgd = ax.legend(handles=[std_x1_patch, std_x2_patch],
                    loc='upper right', fontsize=8,
                    frameon=False)

    # Annotations
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.set_xlabel('Alternative Allele Depth', fontsize=16)
    #ax.set_ylabel('Relative Frequency', fontsize=16)
    ax.set_ylabel('Number of Variants', fontsize=16)
    ax.set_title(os.path.basename(VCF_FILE))

    fig_fname = os.path.basename(VCF_FILE).replace('_report.vcf.gz',
                                                   '_snp_dist.pdf')
    plt.savefig(fig_fname, format='pdf', bbox_inches='tight')
    return


def main(VCF_FILE: str, ASMBL_METADATA: str):
    """
        Plot distribution of depths to work out
        which min_read_number would yield 0-2 snp
    """
    SNP_DEPTH = get_SNP_data(VCF_FILE)
    ASMBL_DEPTH, ASMBL_DEPTH_STD = get_asmbl_depth(ASMBL_METADATA)
    plot_histogram(VCF_FILE, SNP_DEPTH,
                   ASMBL_DEPTH, ASMBL_DEPTH_STD)
    return


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
