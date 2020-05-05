#!/usr/bin/env python3


import toytree
import toyplot
import toyplot.svg
import numpy as np
import argparse

def draw_phylogenetic_tree(filename, root):
    '''
    Takes a phylogenetic tree in Newick format and draw it using
    ToyTree API.
    '''
    tre = toytree.tree(filename, tree_format=0)
    rtre = tre.root(root)
    canvas, axes = rtre.draw(width=600,
                             height=400,
                             node_labels=rtre.get_node_values('support', 0, 0),
                             node_sizes=16,
                             tip_labels_align=True,
                             scalebar=True)
    toyplot.svg.render(canvas, 'tree_plot.svg')

def argument_parser():
    '''Command line argument parser.'''
    parser = argparse.ArgumentParser()
    parser.add_argument('treefile',
                        help='Tree file in Newick format')
    parser.add_argument('root',
                        help='Node used for rooting')
    args = parser.parse_args()
    return args.treefile, args.root

def main():
    tree, root = argument_parser()
    draw_phylogenetic_tree(tree, root)

if __name__ == '__main__':
    main()
