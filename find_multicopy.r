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

            if (x[1,7] == "+") { # due to previous filter steps, all hits are assumed to be on the same strand
                # gap ignorant:
#                 sstart1 <- min(x[i,5:6]) + common_qrange[1] - x[i,2] # for hits on the '-'-strand, BLAST yields reversed sranges!
#                 sstart2 <- min(x[j,5:6]) + common_qrange[1] - x[j,2]
                
                # gap sensitive:
                sstart1 <- system(
                    paste("gawk -v query_ind=", common_qrange[1], " -v query_start=", x[i,2], " -v ref_start=", min(x[i,5:6]), " -f utils/find_pos_btop_plus.awk", sep=""),
                    input = as.character(x[i,16]),
                    intern = T
                )
                sstart2 <- system(
                    paste("gawk -v query_ind=", common_qrange[1], " -v query_start=", x[j,2], " -v ref_start=", min(x[j,5:6]), " -f utils/find_pos_btop_plus.awk", sep=""),
                    input = as.character(x[j,16]),
                    intern = T
                )
            }
            else {
                # gap ignorant:
#                 sstart1 <- min(x[i,5:6]) - common_qrange[1] + x[i,3]
#                 sstart2 <- min(x[j,5:6]) - common_qrange[1] + x[j,3]
                
                # gap sensitive:
                sstart1 <- system(
                    paste("gawk -v query_ind=", common_qrange[1], " -v query_start=", x[i,2], " -v ref_end=", max(x[i,5:6]), " -f utils/find_pos_btop_minus.awk", sep=""),
                    input = as.character(x[i,16]),
                    intern = T
                )
                sstart2 <- system(
                    paste("gawk -v query_ind=", common_qrange[1], " -v query_start=", x[j,2], " -v ref_end=", max(x[j,5:6]), " -f utils/find_pos_btop_minus.awk", sep=""),
                    input = as.character(x[j,16]),
                    intern = T
                )
            }

            sstart1 = as.numeric(sstart1)
            sstart2 = as.numeric(sstart2)
            
            if (sstart1 == -1 || sstart2 == -1) {
                stop(paste("Error: no valid reference position found while processing ", args[1], ", filter_multicopy.r aborted", sep=""))
            }
            
            if (sstart1 != sstart2) {
                mc_length <- common_qrange[2]-common_qrange[1]+1

                if (mc_length > max_mc_length) {
                    cat("removing ", args[1], ": multicopy region found at reference positions ", sstart1, " and ", sstart2, ", length = ", mc_length, " bp\n", sep="", file=args[3], append=TRUE)
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
