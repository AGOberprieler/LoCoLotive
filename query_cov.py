#!/usr/bin/python3
"""
visualize mapped parts of query sequence
usage: ./query_cov.py <marker sequence as single line FASTA> <gff3 file containing query intervals>
"""

import sys

with open(sys.argv[1], "r") as f:
    seq = f.readlines()

with open(sys.argv[2], "r") as f:
    blast_hits = f.readlines()

n = len(seq[1].rstrip())
cov = []

blast_hits = [line.split("\t") for line in blast_hits]

sys.stdout.writelines(seq)

for i in range(len(blast_hits)):
    sys.stdout.write(">" + "-".join(blast_hits[i][3:5]) +"\n")
    sys.stdout.write("".join([
        "-" * (int(blast_hits[i][3]) - 1),
        "*" * (int(blast_hits[i][4]) - int(blast_hits[i][3]) + 1),
        "-" * (n - int(blast_hits[i][4])), "\n"
    ]))
