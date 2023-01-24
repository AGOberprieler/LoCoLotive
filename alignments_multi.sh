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

rm -rf "$outdir"/{genomic_sequences,query_sequences,alignments,query_coverage,temp.fasta,temp.bed,introns.tmp,summary.txt,mafft.log}
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

    # 2) query (possibly reverse-complemented)
    ##########################################

    bedtools getfasta \
        -s \
        -fi "$probes" \
        -bed "${outdir}/query_intervals/${id}.gff3" \
        -fo "${outdir}/query_sequences/${id}.fasta"

    fasta_formatter -i "${outdir}/query_sequences/${id}.fasta" \
        | awk '{ if ($0 !~ /^>/) {print tolower($0)} else {print $0} }' > "${outdir}/temp.fasta"

    mv "${outdir}/temp.fasta" "${outdir}/query_sequences/${id}.fasta"

    echo -e "probe ID: $id\n" >> "${outdir}/mafft.log"

    # align genomic and query sequences
    ###################################
    mafft \
        --inputorder \
        --preservecase \
        --addfragments "${outdir}/query_sequences/${id}.fasta" \
        "${outdir}/genomic_sequences/${id}.fasta" \
    > "${outdir}/alignments/${id}.fasta" 2>> "${outdir}/mafft.log"

    # Sometimes mafft --addframgents fails if the reference contains too much ambiguous nucleotides.
    # --maxambiguous 1.0 could help, but is only available in MAFFT's online version.
    # As a workaround, missing alignments are recomputed without the --addfragments option:
    if ! test -s "${outdir}/alignments/${id}.fasta"; then
        echo -e "\nWarning: mafft --addfragments failed, alignment repeated without the --addfragments option\n" >> "${outdir}/mafft.log"
        cat "${outdir}/genomic_sequences/${id}.fasta" "${outdir}/query_sequences/${id}.fasta" > "${outdir}/mafft_input.tmp"

        mafft \
            --inputorder \
            --preservecase \
            "${outdir}/mafft_input.tmp" \
        > "${outdir}/alignments/${id}.fasta" 2>> "${outdir}/mafft.log"

    fi

    echo -e "\n\n-------------\n\n" >> "${outdir}/mafft.log"

    # visualize query (marker) coverage
    fasta_formatter -i "$probes" | grep "^>${id}$" -A1 > "${outdir}/temp.fasta"
    ./query_cov.py "${outdir}/temp.fasta" "${outdir}/query_intervals/${id}.gff3" > "${outdir}/query_coverage/${id}.fasta"

    progress_bar "$iter" "$i_max" 40 processed
    (( iter++ ))
done


# check for overlapping loci
############################

./detect_overlap.py "${outdir}/genomic_ranges" "${outdir}"

# calculate MSAs for each group of overlapping loci
###################################################

if grep -q "," "${outdir}/groups_of_overlapping_loci.txt"
then
    echo "overlapping loci found, compute group-wise alignments..."

    rm -rf "$outdir/"{genomic_ranges_groupwise,genomic_sequences_groupwise,alignments_groupwise}
    mkdir -p "$outdir/"{genomic_ranges_groupwise,genomic_sequences_groupwise,alignments_groupwise}

    echo -e "\n\n-------------\n-------------\ngroup-wise alignments:\n\n" >> "${outdir}/mafft.log"

    i_max=$(grep -c "^group" "${outdir}/groups_of_overlapping_loci.txt")

    for i_group in $(seq 1 "$i_max"); do
        ids_concatenated=$(sed -n "${i_group}p" "${outdir}/groups_of_overlapping_loci.txt" | sed 's/^group [0-9]\+: //')
        IFS=',' read -r -a ids <<< "$ids_concatenated"

        if [ ${#ids[@]} -eq 1 ]; then

            cp "${outdir}/alignments/${ids[0]}.fasta" "${outdir}/alignments_groupwise/group${i_group}.fasta"

        else

            # extract/modify sequences

            # 1) genome ("+" strand)
            ########################

            for id in "${ids[@]}"; do
                cat "${outdir}/genomic_ranges/${id}.gff3"
            done > "${outdir}/ranges.tmp"

            ref_start=$(cut -f4 "${outdir}/ranges.tmp" | sort -n | head -n1)
            ref_end=$(cut -f5 "${outdir}/ranges.tmp" | sort -n | tail -n1)

            cat "${outdir}/genomic_ranges/${ids[0]}.gff3" | awk -v "i=$ref_start" -v "j=$ref_end" -v quid="$ids_concatenated" '
            BEGIN {
                OFS = "\t"
            }
            {
                $4 = i
                $5 = j
                $9 = quid
                print
            }
            ' > "$outdir/genomic_ranges_groupwise/group${i_group}.gff3"

            bedtools getfasta \
                -fi "$genome" \
                -bed "${outdir}/genomic_ranges_groupwise/group${i_group}.gff3" \
                -fo "${outdir}/genomic_sequences_groupwise/${i_group}.fasta"

            # switch to lower case
            fasta_formatter -i "${outdir}/genomic_sequences_groupwise/${i_group}.fasta" \
                | awk '{ if ($0 !~ />/) {print tolower($0)} else {print $0} }' \
                > "${outdir}/temp.fasta"

            mv "${outdir}/temp.fasta" "${outdir}/genomic_sequences_groupwise/${i_group}.fasta"

            # find and highlight intronic regions
            if [ $# -eq 4 ]
            then
                gff2bed < "$outdir/genomic_ranges_groupwise/group${i_group}.gff3" | cut -f1-6 > "${outdir}/temp.bed"
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

                    cat "${outdir}/genomic_sequences_groupwise/${i_group}.fasta" | utils/upint.py $a $b > "${outdir}/temp.fasta"
                    mv "${outdir}/temp.fasta" "${outdir}/genomic_sequences_groupwise/${i_group}.fasta"
                done
            fi


            echo -e "group ${i_group} (note: group indices can be changed later)\n" >> "${outdir}/mafft.log"

            # 2) combine queries (already extracted)
            ########################################
            for id in "${ids[@]}"; do
                cat "${outdir}/query_sequences/${id}.fasta"
            done > "${outdir}/queries.tmp"

            # align genomic and query sequences
            ###################################
            mafft \
                --inputorder \
                --preservecase \
                --addfragments "${outdir}/queries.tmp" \
                "${outdir}/genomic_sequences_groupwise/${i_group}.fasta" \
            > "${outdir}/alignments_groupwise/group${i_group}.fasta" 2>> "${outdir}/mafft.log"

            # Sometimes mafft --addframgents fails if the reference contains too much ambiguous nucleotides.
            # --maxambiguous 1.0 could help, but is only available in MAFFT's online version.
            # As a workaround, missing alignments are recomputed without the --addfragments option:
            if ! test -s "${outdir}/alignments_groupwise/group${i_group}.fasta"; then
                echo -e "\nWarning: mafft --addfragments failed, alignment repeated without the --addfragments option\n" >> "${outdir}/mafft.log"
                cat "${outdir}/genomic_sequences_groupwise/${i_group}.fasta" "${outdir}/queries.tmp" > "${outdir}/mafft_input.tmp"

                mafft \
                    --inputorder \
                    --preservecase \
                    "${outdir}/mafft_input.tmp" \
                > "${outdir}/alignments_groupwise/group${i_group}.fasta" 2>> "${outdir}/mafft.log"

            fi

            echo -e "\n\n-------------\n\n" >> "${outdir}/mafft.log"

        fi

        progress_bar "$i_group" "$i_max" 40 processed
    done

    rm -rf "$outdir"/{genomic_sequences_groupwise,genomic_ranges_groupwise,ranges.tmp,queries.tmp}
fi


# summarize results
###################

echo -e "\nsummarize results..."
Rscript --vanilla summarize.r "${outdir}"

# add group column (if overlapping loci present)
if grep -q "," "${outdir}/groups_of_overlapping_loci.txt"
then
    cat "${outdir}"/groups_of_overlapping_loci.txt | awk -F "[ ,]" '{gsub(/:/, "", $2); for (i=3; i<=NF; i++) {print $i"\t"$2}}' | sort -k 1b,1 > group_mapping.txt
    sort -k 1b,1 "${outdir}/summary.txt" > "${outdir}/summary.tmp"
    join -t$'\t' -j 1 "${outdir}/summary.tmp" group_mapping.txt > "${outdir}/summary.txt"
fi

# sort table
sort -k3,3nr -k2,2n -k1,1 "${outdir}/summary.txt" > "${outdir}/summary.tmp"
mv "${outdir}/summary.tmp" "${outdir}/summary.txt"

# change group indices for better clarity (according to order of appearance in summary.txt):
if grep -q "," "${outdir}/groups_of_overlapping_loci.txt"
then
    # 1) in tabular summary file
    gawk -v "outdir=$outdir" '
    BEGIN {
        PROCINFO["sorted_in"] = "@ind_num_asc"
        max = 1
        OFS = "\t"
    }
    {
        if (!new[$6]) {
            new[$6] = max
            $6 = max
            max++
        }
        else {
            $6 = new[$6]
        }
        print
    }
    END {
        for (i in new) { 
            printf new[i]" " > outdir"/group_order.tmp"
        }
    }' "$outdir"/summary.txt > "$outdir"/summary.tmp
    mv "$outdir"/summary.tmp "$outdir"/summary.txt

    # 2) in alignment file names
    arr=( $(cat "$outdir"/group_order.tmp) )
    for i in "${!arr[@]}"; do
        old=$(( i + 1 ))
        new="${arr[$i]}"
        mv "$outdir"/alignments_groupwise/group${old}.fasta "$outdir"/alignments_groupwise/group${new}.tmp
    done
    for i in "${!arr[@]}"; do
        old=$(( i + 1 ))
        new="${arr[$i]}"
        mv "$outdir"/alignments_groupwise/group${new}.tmp "$outdir"/alignments_groupwise/group${new}.fasta
    done

    # 3) in list of groups
    awk '{for (i=1; i<=NF; i++) {print "s/^group "i":/group "$i"/"}}' "$outdir"/group_order.tmp > "$outdir"/replacement.sed
    sed -f "$outdir"/replacement.sed "$outdir"/groups_of_overlapping_loci.txt | sort -k2,2n | awk '{$2=$2":"; print}' > "$outdir"/groups_of_overlapping_loci.tmp
    mv "$outdir"/groups_of_overlapping_loci.tmp "$outdir"/groups_of_overlapping_loci.txt

    rm "$outdir"/{group_order.tmp,replacement.sed}
fi

rm -rf "$outdir"/{genomic_sequences,query_sequences,temp.fasta,temp.bed,introns.tmp,mafft_input.tmp}
