#!/bin/bash
usage="
usage: ./filter_angiosperms353.sh <column index> <filter>

possible column indices:
    2 = order
    3 = family
    4 = genus
    5 = species
for possible filter criteria (case-sensitive), see abbreviations.txt

example: ./filter_angiosperms353.sh 2 Fabales > fabales.fasta
"

col="$1"
filter="$2"

if [ $# -ne 2 ]; then
    echo "Error: wrong numer of arguments"
    echo "$usage"
    exit
fi

if [ "$col" -gt 5 ] || [ "$col" -lt 1 ]; then
    echo "Error: column $col does not exist"
    exit
fi

gawk -v "c=$col" -v "f=$filter" '
ARGIND==1 {
    dict[$1] = $c
}
ARGIND==2 {
    if ($0 ~ /^>/) {
        ID = $0
        gsub(/^>/, "", ID)
        gsub(/-.*/, "", ID)
        if (dict[ID] == f) {
            p = 1
        }
        else {
            p = 0
        }
    }
    if (p == 1) {
        print
    }
}
' abbreviations.txt Angiosperms353_targetSequences.fasta
