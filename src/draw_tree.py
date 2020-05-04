#!/usr/bin/env python3

import toytree
import toyplot
import toyplot.png
import numpy as np
import sys

def draw_phylogenetic_tree(filename):
    '''
    Takes a phylogenetic tree in Newick format and draw it using
    ToyTree API.
    '''
    tre = toytree.tree(filename, tree_format=0)
    rtre = tre.root('KF891883.1')
    canvas, axes = rtre.draw(width=600,
                             height=400,
                             node_labels=rtre.get_node_values('support', 0, 0),
                             node_sizes=20,
                             tip_labels_align=True,
                             scalebar=True)
    toyplot.png.render(canvas, 'tree_plot.png')

if __name__ == '__main__':
    draw_phylogenetic_tree(sys.argv[1])
