#!/usr/bin/python3
"""
convert sequence interval to uppercase

input: single line FASTA via stdin
-> run fasta_formatter if necessary
"""

import sys

seq = list(sys.stdin.readlines())

# BED convention: 0-based, half open intervals
i = max(0, int(sys.argv[1]))
j = max(0, int(sys.argv[2]))

if i>j:
    sys.exit("error: start > end")

print("".join((seq[0], seq[1][:i], seq[1][i:j].upper(), seq[1][j:])))
