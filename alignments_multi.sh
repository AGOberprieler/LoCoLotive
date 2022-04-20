#!/bin/bash

usage="usage: ./alignments_multi.sh <probes FASTA> <genome FASTA> <output directory> [<intronic regions BED>]"

if [ $# -lt 3 ]
then
    echo "Error in $0: wrong number of command arguments"
    echo "$usage"
    exit 1
fi

probes=$1
genome=$2
outdir=$3
intronic=$4

source utils/progress.sh

rm -rf "$outdir"/{genomic_sequences,query_sequences,alignments,query_coverage,temp.fasta,temp.bed,introns.tmp,summary.txt}
mkdir -p "$outdir"/{genomic_sequences,query_sequences,alignments,query_coverage}

shopt -s nullglob

iter=1
i_max=$(ls -1q "${outdir}/genomic_ranges" | wc -l)

for f in "${outdir}/genomic_ranges"/*
do
    # extract/modify sequences

    # 1) genome ("+" strand)
    ########################

    id=$(echo "$f" | sed 's/.*\/\([^/]\+\)\.gff3$/\1/')

    bedtools getfasta \
        -fi "$genome" \
        -bed "$f" \
        -fo "${outdir}/genomic_sequences/${id}.fasta"

    # switch to lower case
    fasta_formatter -i "${outdir}/genomic_sequences/${id}.fasta" \
        | awk '{ if ($0 !~ />/) {print tolower($0)} else {print $0} }' \
        > "${outdir}/temp.fasta"

    mv "${outdir}/temp.fasta" "${outdir}/genomic_sequences/${id}.fasta"

    # find and highlight intronic regions
    if [ $# -eq 4 ]
    then
        gff2bed < "$f" | cut -f1-6 > "${outdir}/temp.bed"
        bedtools intersect \
            -a "${outdir}/temp.bed" \
            -b "$intronic" \
            -wo \
            | awk 'BEGIN{OFS="\t"}{print $4,$5,$6,$7,$8,$9}' \
        > "${outdir}/introns.tmp"

        n_introns=$(wc -l < "${outdir}/introns.tmp")

        range_start=$(cut -f2 "${outdir}/temp.bed")
        for k in $(seq 1 "$n_introns")
        do
            read -r intron_start intron_end <<<"$(sed -n "${k}p" "${outdir}/introns.tmp" | cut -f5,6)"
            # intron positions relative to pair
            a=$(( intron_start - range_start ))
            b=$(( intron_end - range_start ))

            cat "${outdir}/genomic_sequences/${id}.fasta" | utils/upint.py $a $b > "${outdir}/temp.fasta"
            mv "${outdir}/temp.fasta" "${outdir}/genomic_sequences/${id}.fasta"
        done
    fi

    # 2) query
    ##########

    bedtools getfasta \
        -s \
        -fi "$probes" \
        -bed "${outdir}/query_intervals/${id}.gff3" \
        -fo "${outdir}/query_sequences/${id}.fasta"

    fasta_formatter -i "${outdir}/query_sequences/${id}.fasta" \
        | awk '{ if ($0 !~ /^>/) {print tolower($0)} else {print $0} }' > "${outdir}/temp.fasta"

    mv "${outdir}/temp.fasta" "${outdir}/query_sequences/${id}.fasta"


    # align genomic and query sequences
    mafft \
        --inputorder \
        --preservecase \
        --addfragments "${outdir}/query_sequences/${id}.fasta" \
        "${outdir}/genomic_sequences/${id}.fasta" \
    > "${outdir}/alignments/${id}.fasta" 2>/dev/null

    # visualize query (marker) coverage
    fasta_formatter -i "$probes" | grep "^>${id}$" -A1 > "${outdir}/temp.fasta"
    ./query_cov.py "${outdir}/temp.fasta" "${outdir}/query_intervals/${id}.gff3" > "${outdir}/query_coverage/${id}.fasta"
    
    progress_bar "$iter" "$i_max" 40 processed
    (( iter++ ))
done


# summarize results
###################

echo -e "\nsummarize results..."
Rscript --vanilla summarize.r "${outdir}"

# sort table
sort -k3 -k2 -nr "${outdir}/summary.txt" > "${outdir}/summary.tmp"
mv "${outdir}/summary.tmp" "${outdir}/summary.txt"

rm -rf "$outdir"/{genomic_sequences,query_sequences,temp.fasta,temp.bed,introns.tmp}
