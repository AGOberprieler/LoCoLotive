#!/usr/bin/env Rscript
library(ape, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE)
library(seqinr, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE)

args = commandArgs(trailingOnly=TRUE)
outdir <- args[1]

# remove old summary file
unlink(file.path(outdir, "summary.txt"), expand = F)

# IDs of target sequences having passed all filtering steps
IDs <- list.files(path=file.path(outdir, "hits_filtered"), all.files=T, full.names=F, no..=T)

for (ID in IDs) {

    # start index (1-based) of genomic fragment covering all BLAST hits for a given target sequence
    genomic_start <- read.csv(
        file.path(outdir, "genomic_ranges", paste(ID, ".gff3", sep="")),
        header = F,
        sep = "\t"
    )[,4]
    genomic_start <- as.numeric(genomic_start)

    # blast_hits <- read.csv(args[1], sep="\t", header=F)
    
    # read BLAST results
    blast_hits <- read.csv(
        file.path(outdir, "hits_filtered", ID),
        sep = "\t",
        header = F
    )

    # reformat interval indices following GFF3 conventions
    gff_intervals <- t(apply(as.matrix(blast_hits[,5:6]), 1, sort))
    colnames(gff_intervals) <- NULL
    #strands <- strands[order(gff_intervals[,1]), ]
    gff_intervals <- gff_intervals[order(gff_intervals[,1]), ]

    # ref_seq <- args[2]
    # ref_offset <- as.numeric(args[3])

    # seq <- read.fasta(ref_seq, seqtype = "DNA", as.string = FALSE, forceDNAtolower = FALSE)[[1]]
    
    # read genomic sequence (where intronic regions are possibly highlighted as uppercase)
    genomic_seq <- read.fasta(
        file.path(outdir, "genomic_sequences", paste(ID, ".fasta", sep="")),
        seqtype = "DNA",
        as.string = FALSE,
        forceDNAtolower = FALSE
    )[[1]]

    # count number of (intronic) base pairs between consecutive BLAST hits
    bp_between <- NULL
    intronic_between <- NULL
    
    for (i in seq_along(gff_intervals[,1])) {
        if (i > 1) {
            curr_bp <- gff_intervals[i,1] - end - 1
            if (curr_bp > 0) {
                curr_intronic <- sum(grepl("[A-Z]", genomic_seq[(end + 2 - genomic_start) : (gff_intervals[i,1] - genomic_start)]))
            }
            else {
                curr_intronic <- 0
            }
            bp_between <- c(bp_between, curr_bp)
            intronic_between <- c(intronic_between, curr_intronic)
        }
        end <- gff_intervals[i,2]
    }
    
    # alignment <- read.dna(args[4], format = "fasta", as.character = TRUE, as.matrix = TRUE)
    
    # get alignment length
    alignment <- read.dna(
        file.path(outdir, "alignments", paste(ID, ".fasta", sep="")),
        format = "fasta",
        as.character = TRUE,
        as.matrix = FALSE
    )

    ali_length <- length(alignment[[1]])
    n_hits <- length(alignment) - 1
    # ID <- sub(".*/", "", args[1])

    # append results to summary file
    cat(
        ID,
        ali_length,
        n_hits,
        paste(bp_between, collapse=","),
        paste(intronic_between, collapse=","),
        sep = "\t",
        file = file.path(outdir, "summary.txt"),
        append = T
    )
    cat("\n", file = file.path(outdir, "summary.txt"), append = T)
}
