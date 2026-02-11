#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import sys
import os


def extract_depth(FNAME: str) -> list:
    """
        Read FNAME and return a numpy array with
        depth data per nucleotide of the reference genome.
    """
    data = list()
    for entry in open(FNAME, 'r'):
        entry = entry.strip('\n')
        # Extract data
        reference, pos, depth = entry.split('\t')

        # Append to list
        data.append(int(depth))
    return np.array(data)


def histogram_plot(HIST_DATA: list, FNAME: str):
    """
        Plot histogram
    """
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6,5))

    # Plot data
    ax.hist(HIST_DATA, bins=250, width=0.5, color='black')
    y_lims = ax.get_ylim()

    # Plot accessory data
    mean = HIST_DATA.mean()
    std = HIST_DATA.std()
    ax.plot([mean, mean], y_lims, color='red',
            zorder=-1, linewidth=4)
    std_x1_patch = plt.Rectangle(xy=[mean - std, y_lims[0]],
                                 width=std*2,  # ± std
                                 height=y_lims[1],
                                 color='grey', alpha=0.5,
                                 zorder=-5,
                                 label='Mean ± 1 x StDev')
    std_x2_patch = plt.Rectangle(xy=[mean - (std*2), y_lims[0]],
                                 width=std*2*2,  # ± 2 x std
                                 height=y_lims[1],
                                 color='lightgrey', alpha=0.5,
                                 zorder=-5,
                                 label='Mean ± 2 x StDev')
    ax.add_patch(std_x1_patch)
    ax.add_patch(std_x2_patch)

    # Lebend
    lgd = ax.legend(handles=[std_x1_patch, std_x2_patch], loc='upper right',
                    fontsize=8, frameon=False, bbox_to_anchor=(1, 1.075),
                    ncols=2)

    # Labels
    ax.set_xlabel('Depth', fontsize=18)
    ax.set_ylabel('Number of basepairs', fontsize=18)
    ax.tick_params(axis='both', labelsize=14)
    
    # Format axes
    ax.set_ylim(y_lims)  # Restore original limit, since adding data will change it
    ax.ticklabel_format(axis='y', style='scientific',
                        scilimits=(0,0), useMathText=True)

    # Save
    save_fig = os.path.basename(FNAME).replace('_depth.tsv', '_depth_hist.pdf')
    plt.savefig(save_fig, format='pdf', bbox_inches='tight')
    return


def main(FNAME):
    depth_data = extract_depth(FNAME)
    histogram_plot(depth_data, FNAME)



if __name__ == '__main__':
    main(sys.argv[1])
