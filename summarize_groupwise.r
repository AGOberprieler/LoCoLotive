#!/usr/bin/env Rscript
library(ape, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE)
library(seqinr, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE)

args = commandArgs(trailingOnly=TRUE)
outdir <- args[1]

# remove old summary file
unlink(file.path(outdir, "summary_groupwise.txt"), expand = F)

groups <- strsplit(
    read.csv(file.path(outdir, "groups_of_overlapping_loci.txt"), sep=" ", header=F, stringsAsFactors=F)[,3],
    ","
)

for (group_ID in seq_along(groups)) {

    if (length(groups[[group_ID]]) == 1) {
        # start index (1-based) of genomic fragment covering all BLAST hits for a given probe sequence
        genomic_start <- read.csv(
            file.path(outdir, "genomic_ranges", paste(groups[[group_ID]], ".gff3", sep="")),
            header = F,
            sep = "\t"
        )[,4]
    }
    else {
        # start index (1-based) of genomic fragment covering all BLAST hits for a given group of probe sequences
        genomic_start <- read.csv(
            file.path(outdir, "genomic_ranges_groupwise", paste("group", group_ID, ".gff3", sep="")),
            header = F,
            sep = "\t"
        )[,4]
    }
    genomic_start <- as.numeric(genomic_start)

    # read BLAST results
    blast_hits <- NULL
    for (probe_ID in groups[[group_ID]]) {
        blast_hits_probe <- read.csv(
            file.path(outdir, "hits_filtered", probe_ID),
            sep = "\t",
            header = F
        )
        blast_hits <- rbind(blast_hits, blast_hits_probe)
    }

    # reformat interval indices following GFF3 conventions
    gff_intervals <- t(apply(as.matrix(blast_hits[,5:6]), 1, sort))
    colnames(gff_intervals) <- NULL
    #strands <- strands[order(gff_intervals[,1]), ]
    gff_intervals <- gff_intervals[order(gff_intervals[,1]), ]

    # merge overlapping intervals
    intervalls_merged <- NULL

    for (i in seq_along(gff_intervals[,1])) {
        if (i==1) {
            start <- gff_intervals[i,1]
            end <- gff_intervals[i,2]
        }
        else {
            # overlapping
            if (gff_intervals[i,1] <= end) {
                end <- max(end, gff_intervals[i,2])
            }
            # non-overlapping
            else {
                intervalls_merged <- rbind(intervalls_merged, c(start, end))
                start <- gff_intervals[i,1]
                end <- gff_intervals[i,2]
            }
        }
    }

    intervalls_merged <- rbind(intervalls_merged, c(start, end))
    gff_intervals <- intervalls_merged

    # read genomic sequence (where intronic regions are possibly highlighted as uppercase)
    if (length(groups[[group_ID]]) == 1) {
        genomic_seq <- read.fasta(
            file.path(outdir, "genomic_sequences", paste(groups[[group_ID]], ".fasta", sep="")),
            seqtype = "DNA",
            as.string = FALSE,
            forceDNAtolower = FALSE
        )[[1]]
    }
    else {
        genomic_seq <- read.fasta(
            file.path(outdir, "genomic_sequences_groupwise", paste(group_ID, ".fasta", sep="")),
            seqtype = "DNA",
            as.string = FALSE,
            forceDNAtolower = FALSE
        )[[1]]
    }

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

    # get alignment length
    alignment <- read.dna(
        file.path(outdir, "alignments_groupwise", paste("group", group_ID, ".fasta", sep="")),
        format = "fasta",
        as.character = TRUE,
        as.matrix = FALSE
    )

    ali_length <- length(alignment[[1]])
    # n_hits <- length(alignment) - 1
    n_hit_groups <- nrow(gff_intervals)

    # append results to summary file
    cat(
        group_ID,
        ali_length,
        n_hit_groups,
        paste(bp_between, collapse=","),
        paste(intronic_between, collapse=","),
        sep = "\t",
        file = file.path(outdir, "summary_groupwise.txt"),
        append = T
    )
    cat("\n", file = file.path(outdir, "summary_groupwise.txt"), append = T)
}
