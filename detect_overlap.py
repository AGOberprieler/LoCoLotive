#!/usr/bin/env python3
"""
This script detects overlaps between loci by looking at their location on the reference genome
and then identifies groups of overlapping loci.
Usage: ./detect_overlap.py <input directory> <output directory>
The input directory has to comprise genomic ranges in gff3 format (one file per locus, informative file names).
"""

from collections import namedtuple
import glob
import os
import sys
import numpy as np
from scipy.sparse.csgraph import connected_components

region = namedtuple("region", ["reference", "start", "end"])

if not os.path.isdir(sys.argv[1]):
    raise ValueError(f'Error: directory "{sys.argv[1]}" does not exit!')
if not os.path.isdir(sys.argv[2]):
    raise ValueError(f'Error: directory "{sys.argv[2]}" does not exit!')


dict_IDs = dict()
overlapping = dict()

# read input files
for path in glob.glob(os.path.join(sys.argv[1], "*gff3")):
    fname = os.path.basename(path)
    ID = os.path.splitext(fname)[0]

    with open(path, "r") as infile:
        cols = infile.readline().split("\t")
        if len(cols) != 9:
            raise ValueError(f"Error: can't interprete {path}, wrong number of tab-delimited columns!")

        dict_IDs[ID] = region(cols[0], cols[3], cols[4])
        overlapping[ID] = []

adj_mat = np.zeros((len(dict_IDs), len(dict_IDs)), dtype=np.bool_)
IDs = sorted(dict_IDs.keys())

# identify pairs of overlapping loci
for i in range(len(IDs)):
    for j in range(i+1, len(IDs)):
        id1 = IDs[i]
        id2 = IDs[j]
        region1 = dict_IDs[id1]
        region2 = dict_IDs[id2]

        if region1.reference == region2.reference:
            if (
                (region1.start <= region2.end and region1.start >= region2.start) or
                (region2.start <= region1.end and region2.start >= region1.start)
            ):
                overlapping[id1].append(id2)
                overlapping[id2].append(id1)
                adj_mat[i,j] = 1
                adj_mat[j,i] = 1


# report overlapping loci (locus-wise)
with open(os.path.join(sys.argv[2], "overlapping_loci.txt"), "w") as outfile:
    for ID in IDs:
        outfile.write(ID + ": " + ",".join(overlapping[ID]) + "\n")


# identify groups of overlapping loci
n_comp, comp_labels = connected_components(adj_mat, directed=False)

# report groups
with open(os.path.join(sys.argv[2], "groups_of_overlapping_loci.txt"), "w") as outfile:
    for i_comp in range(n_comp):
        inds_group = np.where(comp_labels == i_comp)[0].tolist()
        outfile.write(f"group {i_comp + 1}: " + ",".join([IDs[i] for i in inds_group]) + "\n")
