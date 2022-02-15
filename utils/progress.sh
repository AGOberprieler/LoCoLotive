#!/bin/bash

# Print simple progress bar.
# usage: progress_bar <current iteration> <total number of iterations> <bar width> <unit>
# The first three arguments are required, <unit> may be missing.
progress_bar() {
    nchar1=$(( $3 * $1 / $2 ))
    nchar2=$(( $3 - nchar1 ))

    printf "\r["
    printf "%${nchar1}s" | tr " " "#"
    printf "%${nchar2}s] $1 of $2 $4"
    
    if [ $nchar2 -eq 0 ]; then
        printf "\n"
    fi
}
