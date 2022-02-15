#!/bin/bash

usage="usage: ./extract_introns.sh <annotation GFF3> <output directory>"

if [ $# -ne 2 ]
then
    echo "Error in $0: wrong number of command arguments"
    echo "$usage"
    exit 1
fi

annotation=$1
outdir=$2

gt gff3 -addintrons "$annotation" | grep -v "^#" | awk '$3=="intron"' | sort -k1,1 -k4,4n | bedtools merge > "$outdir"/intronic.BED
grep -v "^#" "$annotation" | awk '$3=="exon"' | sort -k1,1 -k4,4n | bedtools merge > "$outdir"/exonic.BED

bedtools subtract -a "$outdir"/intronic.BED -b "$outdir"/exonic.BED > "$outdir"/intronic_strict.BED
