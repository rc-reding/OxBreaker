#!/usr/bin/env python3

import os
import sys
import subprocess
import numpy as np
from ete3 import NodeStyle, Tree, TreeStyle, TextFace
import matplotlib.pyplot as plt


def _get_sequence_type(PATH, barcode) -> str:
    """
        Retrieve sequence type given a sample barcode
    """
    print(PATH, barcode)
    barcode_fname = [f for f in os.listdir(PATH) if barcode.split('.')[0] == f.split('.')[0].rstrip('_ST')]
    print("Barcode fname:", barcode_fname)
    assert len(barcode_fname) == 1

    CMD = str('cat ') + PATH + barcode_fname[0] + str(' | cut -f 3')
    output = subprocess.check_output(CMD, shell=True).decode()
    return output.replace('\n', '')


def render_tree(tree: Tree, MLST_PATH: str, cladogram: bool = False, DEBUG: bool = False) -> Tree:
    """
        Configure visual options for the tree
    """
    # Re-name leafs
    for leaf in tree.iter_leaves():
        seq_type = _get_sequence_type(MLST_PATH, leaf.name)
        leaf.name = leaf.name.split('.')[0].replace('barcode', '').capitalize() + \
            str(' (') + seq_type + str(')')
    # Re-format node style
    for node in tree.traverse():
        node.img_style['size'] = 0  # Hide node circles
        if node.is_leaf():
            node.img_style['size'] = 3
            name_face = TextFace(node.name, fgcolor='k', fsize=10)
            name_face.margin_left = 2.5
            node.add_face(name_face, column=0, position='branch-right')

    ts = TreeStyle()
    ns = NodeStyle()

    # Tree style
    ts.show_leaf_name = False
    ts.margin_bottom = 10
    ts.margin_top = 10
    ts.margin_left = 10
    ts.margin_right = 10

    if DEBUG is False:
        # Default to circular, unrooted trees
        # unless reference is included
        ts.mode = 'c'

    if cladogram is True:
        ts.show_scale = False
        ts.force_topology = True
        ts.show_branch_support = True
        ns.size = 1
    else:
        ts.show_branch_support = False
        #ts.scale_length = 100

    # Node style
    ns['vt_line_width'] = 0.5
    ns['hz_line_width'] = 0.5
    return tree, ts


def fix_notation(TREE_FNAME: str, NDIGITS: int = 4):
    """
    """
    broken_tree = open(TREE_FNAME, 'r').read()
    NDIGITS = str(NDIGITS)
    if broken_tree.count('[') == 0:
        print("Notation correct, nothing to fix.\n")
        return
    else:
        fixed_tree = str('')
        for brackets in range(0, broken_tree.count('[')):
            # Isolate SH-like support
            tmp_tree, sh_val = broken_tree.split(']')[brackets].split('[')
            # Break tmp_tree by ':', then insert sh_val
            tmp_tree = tmp_tree.split(':')
            tmp_tree.insert(-1, sh_val)
            final_fix = list()
            for i in tmp_tree:
                if i.count('.') > 0:
                    if i.count(',') == 1:
                        val, rest = i.split(',')
                        val = str(float(val)).format(':' + NDIGITS + 'e')
                        i = val + ',' + rest
                    else:
                        if len(i.split(')')) > 1:
                            val, rest = i.split(')')
                            val = str(float(val)).format(':' + NDIGITS + 'e')
                            i = val + ')' + rest
                        else:
                            if str(float(i)).format(':' + NDIGITS + 'e') == i:
                                i = str(round(float(i), int(NDIGITS)))
                            else:
                                i = str(float(i)).format(':' + NDIGITS + 'e')
                final_fix.append(i)
            # Rejoin and grow fixed tree
            fixed_tree += str(':').join(final_fix)

        # Add existing data past the latest bracket
        fixed_tree += broken_tree[broken_tree.rfind(']') + 1:]

        # Remove colon after parenthesis:
        fixed_tree = fixed_tree.replace('):', ')')

        # Overwrite tree:
        with open(TREE_FNAME, 'w') as fOut:
            fOut.writelines(fixed_tree)
            fOut.close()
    return


def load_phylogeny(PATH: str, FNAME: str) -> Tree:
    """

    """
    # Fix tree notation before loading, if needed
    # fix_notation(PATH + FNAME, NDIGITS=4)

    # Load tree
    tree = Tree(PATH + FNAME, quoted_node_names=True, format=1)
    DEBUG = False
    if str('debug') in PATH.lower():
        tree.set_outgroup('reference')
        DEBUG = True
    tree.ladderize(direction=1)  # direction=0 means outgroup on TOP
    TREE_ORDER = [entry.split('.')[0].replace('barcode', '').capitalize() for entry in tree.get_leaf_names()]
    return tree, TREE_ORDER, DEBUG


def _rearrange_distance_matrix(DISTANCES: list, LABELS: list, ORDER: list) -> list:
    """
       Re-arrange distance matrix to follow the order in phylogeny for
       downstream plotting.
    """
    TMP_D = list()
    TMP_LAB = list()
    for entry in ORDER:
        try:
            IDX = LABELS.index(entry)
            # Append to new list
            TMP_D.append(DISTANCES[IDX])
            TMP_LAB.append(LABELS[IDX])
        except ValueError:
            pass
    return TMP_D, TMP_LAB


def load_distance_matrix(PATH: str, TREE_ORDER: list, ALL_V_ALL=False) -> np.matrix:
    """
        Load distance matrix from file.
    """
    import csv
    FNAME = str('distance_matrix.tsv')  # All assemblies listed
    CSV_READER = csv.reader(open(PATH + FNAME), delimiter='\t')
    D = list()
    for ROW in CSV_READER:
        if str("snp-dists") in ROW[0] or len(ROW) == 1:
            LABELS = list()
            # for entry in ROW[1:]:
            #     LABELS.append(entry.split('.')[0].replace('barcode',
            #                                               '').capitalize())
        else:
            LABELS.append(ROW[0].split('.')[0].replace('barcode', '').capitalize())
            if ALL_V_ALL is True:
                D.append(ROW[1:])
            else:
                D.append(ROW[1])

    # Re-arrange distances and labels
    if ALL_V_ALL is False:
        D, LABELS = _rearrange_distance_matrix(D, LABELS, TREE_ORDER)

    if ALL_V_ALL is True:
        return np.matrix(D, dtype='int'), LABELS
    else:
        return np.array(D, dtype='int'), LABELS


# Define PATH
PATH = sys.argv[1]
DEST_PATH = sys.argv[2]
MLST_PATH = sys.argv[3]

phylogeny, tree_order, debug_flag = load_phylogeny(PATH, 'corrected_tree.nwk')
Dm, Dm_labels = load_distance_matrix(PATH, tree_order, ALL_V_ALL=False)

tree, ts = render_tree(phylogeny, MLST_PATH, DEBUG=debug_flag)
tree.render(DEST_PATH + 'phylogeny.pdf', tree_style=ts)

## phylogeny, _ = load_phylogeny(PATH, 'sh_tree.nwk')  # Reload to purge previous configuration
#phylogeny, _ = load_phylogeny(PATH, 'corrected_tree.nwk')  # Reload to purge previous configuration
#tree, ts = render_tree(phylogeny, MLST_PATH, cladogram=True)
#tree.render(DEST_PATH + 'cladogram.pdf', tree_style=ts)

