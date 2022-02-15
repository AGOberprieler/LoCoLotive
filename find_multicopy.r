#!/usr/bin/env Rscript
# delete files/markers with at least one multicopy region


args = commandArgs(trailingOnly=TRUE)
x <- read.csv(args[1], header=F, sep="\t")
max_mc_length <- as.numeric(args[2])

# sort by query start
x <- x[order(x[,2]),]

mc_found <- FALSE

for (i in 1:(nrow(x)-1)) {
    for (j in (i+1):nrow(x)) {
        if (x[j,2] <= x[i,3]) {
            common_qrange <- c(x[j,2], min(x[i,3], x[j,3]))

            sstart1 <- min(x[i,5:6]) + common_qrange[1] - x[i,2] # for hits on the '-'-strand, BLAST yields reversed sranges!
            sstart2 <- min(x[j,5:6]) + common_qrange[1] - x[j,2]

            if (sstart1 != sstart2) {
                mc_length <- common_qrange[2]-common_qrange[1]+1

                if (mc_length > max_mc_length) {
                    cat("removing ", args[1], ": multicopy region found, length=", mc_length, "\n", sep="", file=args[3], append=TRUE)
                    mc_found <- TRUE
                    break
                }
            }
        }
    }
    if (mc_found) break
}

if (mc_found) {
    cat("1")
} else {
    cat("0")
}
