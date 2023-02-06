#!/bin/bash
# run preprocess.sh first
# input file: col.1=qid, col.2=qstart, col.3=qend, col.4=chromosome, col.5=sstart, col.6=send, col.7=sstrand(+/-)

usage="usage: ./filter.sh <input file> <output directory> <mc_length>"

if [ $# -ne 3 ]
then
    echo "Error in $0: wrong number of command arguments"
    echo "$usage"
    exit 1
fi

infile=$1
outdir=$2
mc_length=$3

source utils/progress.sh

rm -rf "${outdir}/"{hits_filtered,query_intervals,genomic_ranges,filtering.log}
mkdir -p "${outdir}/"{hits_filtered,query_intervals,genomic_ranges}

# split blast hits by query ID
tail -n+2 "$infile" \
    | awk -v "outdir=${outdir}" '
    BEGIN {
        OFS="\t"
    }
    $7=="plus" {
        $7="+"
    }
    $7=="minus" {
        $7="-"
    }
    {
        print >> outdir"/hits_filtered/"$1
    }'

# Filter target sequences:

shopt -s nullglob

echo "remove query IDs with less than 2 blast hits..."
iter=1
i_max=$(ls -1q "${outdir}/hits_filtered" | wc -l)
n_removed=0
for f in "${outdir}/hits_filtered"/*
do
    # count number of lines
    n=$(wc -l < "$f")
    if [ "$n" -lt 2 ]
    then
        echo "removing $f: less than 2 blast hits" >> "${outdir}/filtering.log"
        rm -f "$f"
        (( n_removed++ ))
    fi

    if [ $(( (i_max - iter) % 25 )) -eq 0 ]
    then
        progress_bar "$iter" "$i_max" 40 processed
    fi
    (( iter++ ))
done
echo -e "$n_removed discarded\n"


echo "remove query IDs with hits on different chromosomes/scaffolds/contigs etc. ..."
iter=1
i_max=$(ls -1q "${outdir}/hits_filtered" | wc -l)
n_removed=0
for f in "${outdir}/hits_filtered"/*
do
    # count number of chromosomes (col. 4)
    n=$(cut -f4 "$f" | sort | uniq | wc -l)
    if [ "$n" -gt 1 ]
    then
        echo "removing $f: blast hits on $n different genomic sequences" >> "${outdir}/filtering.log"
        rm -f "$f"
        (( n_removed++ ))
    fi

    if [ $(( (i_max - iter) % 25 )) -eq 0 ]
    then
        progress_bar "$iter" "$i_max" 40 processed
    fi
    (( iter++ ))
done
echo -e "$n_removed discarded\n"


echo "remove query IDs with hits on both strands..."
iter=1
i_max=$(ls -1q "${outdir}/hits_filtered" | wc -l)
n_removed=0
for f in "${outdir}/hits_filtered"/*
do
    # count number of strands (col. 7)
    n=$(cut -f7 "$f" | sort | uniq | wc -l)
    if [ "$n" -gt 1 ]
    then
        echo "removing $f: blast hits on $n different strands"  >> "${outdir}/filtering.log"
        rm -f "$f"
        (( n_removed++ ))
    fi

    if [ $(( (i_max - iter) % 25 )) -eq 0 ]
    then
        progress_bar "$iter" "$i_max" 40 processed
    fi
    (( iter++ ))
done
echo -e "$n_removed discarded\n"


echo "remove query IDs with multi-copy regions of at least ${mc_length} bp..."
iter=1
i_max=$(ls -1q "${outdir}/hits_filtered" | wc -l)
n_removed=0
for f in "${outdir}/hits_filtered"/*
do
    mc_found=$(Rscript --vanilla find_multicopy.r "$f" "$mc_length" "${outdir}/filtering.log")

    if [ "$mc_found" -eq 1 ]
    then
        rm -f "$f"
        (( n_removed++ ))
    fi

    progress_bar "$iter" "$i_max" 40 processed
     (( iter++ ))
done
echo -e "$n_removed discarded\n"

shopt -u nullglob

if [ -z "$(ls -A "${outdir}"/hits_filtered)" ]; then
    echo "Warning: no BLAST hits left after filtering"
    touch "${outdir}/summary.txt"
    exit
fi

# print info
echo "number of BLAST hits per target sequence after filtering:"
wc -l "${outdir}/hits_filtered"/* | head -n-1 | awk '{print $1}' | sort -n | uniq -c | awk '{print $2" hits: "$1" targets"}'

# summarize filtered blast hits, prepare input for alignments_multi.sh

grep "# Fields:" "$infile" | head -n1  > "${outdir}/hits_filtered.csv" # write header
cat "$outdir"/hits_filtered/* >> "${outdir}/hits_filtered.csv"

shopt -s nullglob
for f in "${outdir}/hits_filtered"/*
do
    # reformat matched query intervals:

    id=${f/*\//}
    awk 'BEGIN{OFS="\t"} {print $1, ".", ".", $2, $3, "0", $7, ".", "."}' "$f" | sort -k4,5n > "${outdir}/query_intervals/${id}.gff3"

    # get genomic range:
    # swap sstart/send if strand = '-'
    awk 'BEGIN {OFS="\t"} {if ($7=="-") {i=$5; $5=$6; $6=i;} print}' "$f" > "${outdir}/range.tmp"
    range_chr=$(cut -f4 "${outdir}/range.tmp" | head -n1)
    range_strand=$(cut -f7 "${outdir}/range.tmp" | head -n1)
    range_start=$(cut -f5 "${outdir}/range.tmp" | sort -n | head -n1)
    range_end=$(cut -f6 "${outdir}/range.tmp" | sort -n | tail -n1)
    echo -e "${range_chr}\t.\t.\t${range_start}\t${range_end}\t0\t${range_strand}\t.\tquid=$id" > "${outdir}/genomic_ranges/${id}.gff3"
done
shopt -u nullglob

rm -f "${outdir}/range.tmp"
