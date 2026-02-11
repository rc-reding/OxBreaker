#!/usr/bin/env python3
import sys
import numpy as np
from matplotlib import pyplot as plt


def load_matrix(MATRIX_FNAME: str):
    """
        Load all vs all SNP distance matrix,
        and return bottom half.
    """
    # Load file. Note first row has number of samples,
    # and last row is empty.
    file_content = open(MATRIX_FNAME, 'r').read().split('\n')[1:-1]
    # matrix_data = dict()
    snp_data = list()
    sample_id = list()
    for row in file_content:
        # sample_id = row.split('\t')[0]
        if len(row.split('\t')) > 1:
            sample_name = row.split('\t')[0].replace('barcode', '')
            sample_name = sample_name.replace('.vcf.gz', '')
            sample_name = sample_name.replace('reference', 'Ref.')
            sample_id.append(sample_name.capitalize())
        data = [int(i) for i in row.split('\t')[1:]]
        # matrix_data[sample_id] = data
        snp_data.append(data)

    # Extract only the triangular/lower section
    # of the distance matrix (includes diagonal)
    idx = np.tril_indices(len(data), -1)
    return np.array(snp_data)[idx], np.matrix(snp_data), tuple(sample_id)


def plot_snp_histogram(matrix_data: np.array):
    """
    """
    fig = plt.figure(figsize=(5, 5))
    plt.hist(matrix_data, bins=50, rwidth=0.95)
    plt.tight_layout()
    
    # Labels
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlabel('No. of SNVs', fontsize=14)
    plt.ylabel('No. of Isolates', fontsize=14)
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0),
                         useMathText=True)
    plt.savefig('snp_histogram.pdf', format='pdf',
                bbox_inches='tight')
    return


def plot_distance_matrix(snp_data: np.array, sample_id: tuple):
    """
    """
    fig = plt.figure(figsize=(5, 5.5))

    # Since matrix is symmetric, extract triangular
    # lower section of the matrix (including diag.)
    idx = np.triu_indices(len(snp_data), 0)
    snp_alpha = np.ones_like(snp_data, dtype='float')
    snp_alpha[idx] = 0.0

    # Plot
    p = plt.pcolor(snp_data, alpha=snp_alpha, edgecolors='black',
                   linewidth=0.1, cmap='OrRd')
    ax = plt.gca()  # Get axis
    
    # Remove frame
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    # Major axis (for labels)
    ax.tick_params(which='major', labeltop=True, length=0,
                   top=True, bottom=False, labelbottom=False)
    plt.xticks(ticks=np.arange(len(sample_id)) + 0.5,
               labels=sample_id, rotation=45,
               fontsize=6, horizontalalignment='left')
    plt.yticks(ticks=np.arange(len(sample_id)) + 0.5,
               labels=sample_id, rotation=0,
               fontsize=6)

    # Annotation
    for (x, y), val in np.ndenumerate(snp_data):
        if snp_alpha[y, x] != 0:
            if val <= 40:
                ax.text(x+0.5, y+0.5, val, fontsize=4,
                        rotation=0,
                        horizontalalignment='center',
                        verticalalignment='center',
                        color='red', weight='bold')
            elif 100 >= val > 40:
                ax.text(x+0.5, y+0.5, val, fontsize=4,
                        rotation=0,
                        horizontalalignment='center',
                        verticalalignment='center',
                        color='black', weight='bold')
            else:
                new_val = str('$>$100')
                ax.text(x+0.5, y+0.5, new_val, fontsize=3,
                        rotation=45,
                        horizontalalignment='center',
                        verticalalignment='center',
                        color='black')

    cb = fig.colorbar(p, ax=ax, extend='max', label='Number of Variants',
                      orientation='horizontal', aspect=50, shrink=0.75, pad=-0.09)
    cb.ax.tick_params(axis='both', direction='out')

    # Labels
    plt.xlabel('Sample ID', fontsize=14)
    plt.ylabel('Sample ID', fontsize=14)
    ax.xaxis.set_label_position('top')

    plt.savefig('distance_matrix.pdf', format='pdf',
                bbox_inches='tight')
    return


def main(FNAME: str):
    """
    
    """
    # Load data
    mx_data, snp_data, sample_id = load_matrix(FNAME)
    print("Matrix loaded.")

    # Plot
    plot_snp_histogram(mx_data)
    print("Histogram done.")
    plot_distance_matrix(snp_data, sample_id)
    print("Matrix done.")
    plt.close('all')

if __name__ == '__main__':
    main(sys.argv[1])
