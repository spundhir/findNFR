#!/usr/bin/env Rscript
#suppressPackageStartupMessages(library("optparse", lib.loc = "/home/users/sachin/R/x86_64-redhat-linux-gnu-library/2.15"))
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
  make_option(c("-i", "--inputDeeptoolsMatrix"), help="input deepPlot matrix file"),
  make_option(c("-j", "--sampleAnno"), help="sample annotation, aveage will be computed for samples with identical name (eg. wt,wt,ko,ko)"),
  make_option(c("-o", "--outputDeeptoolsMatrix"), help="output deepPlot matrix file (should be .gz)")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if((is.null(opt$inputDeeptoolsMatrix) | is.null(opt$outputDeeptoolsMatrix))) {
  cat("\nProgram: deepPlot.R (edit deepPlot matrix file)\n")
  cat("Author: BRIC, University of Copenhagen, Denmark\n")
  cat("Version: 1.0\n")
  cat("Contact: pundhir@binf.ku.dk\n");
  print_help(parser)
  q()
}

suppressPackageStartupMessages(library("jsonlite"))

## read input deepTools matrix file
df <- read.table(pipe(sprintf("zgrep -v '@' %s", opt$inputDeeptoolsMatrix)))
mat <- as.matrix(df[,-c(1:6)])

## read JSON paramters
h <- scan(opt$inputDeeptoolsMatrix, n=1, sep="\n", what=character())
params <- fromJSON(txt=gsub("@", "", h))

## enlist signal matrices
len <- params$sample_boundaries[2]-params$sample_boundaries[1]
mat.lst <- lapply(seq(1,length(params$sample_boundaries)-1), function(x) {
    mat[,c((params$sample_boundaries[x]+1):(params$sample_boundaries[x]+len))]
})

## compute mean signal
if(!is.null(opt$sampleAnno)) {
    samples <- unlist(strsplit(opt$sampleAnno, ","))
    index <- lapply(unique(samples), function(x) grep(x, samples))

    meanMat <- lapply(index, function(x) {
        X <- mat.lst[x]
        Y <- do.call(cbind, X)
        Y <- array(Y, dim=c(dim(mat.lst[[1]]), length(X)))
        colMeans(aperm(Y, c(3, 1, 2)), na.rm = TRUE)
    })
    ## test, if mean works
    # x  <- array(1:24, 2:4)
    # colMeans(aperm(x[,,1:2], c(3,1,2)), na.rm = T)

    meanMat <- do.call(cbind, meanMat)
    #dim(meanMat)

    ## reformat parameters
    meanParams <- lapply(names(params), function(x) {
            if(length(params[[x]]) >= length(mat.lst)) {
                params[[x]][c(1:(length(unique(samples))))]
            } else if(!is.null(params[[x]])) {
                params[[x]]
            } else { "null" }
        })
    names(meanParams) <- names(params)
    meanParams$group_labels <- params$group_labels
    meanParams$group_boundaries <- params$group_boundaries
    meanParams$sample_labels <- unique(samples)
    meanParams$sample_boundaries <- params$sample_boundaries[c(1:(length(meanParams$sample_labels)+1))]
    meanParamsStr <- gsub("\\[\"null\"\\]", "null", toJSON(meanParams))
    #meanParamsStr

    ## some values in the deepTools matrix are "nan"
    meanMat[is.na(meanMat)] <- 0

    gz1 <- gzfile(opt$outputDeeptoolsMatrix, "w")
    sink(gz1)
    cat(sprintf("@%s\n", meanParamsStr))
    write.table(cbind(df[,c(1:6)], meanMat), "", sep="\t", quote = F, row.names = F, col.names = F)
    #write.table(cbind(df[,c(1:6)], mat.lst[[1]], mat.lst[[2]]), "", sep="\t", quote = F, row.names = F, col.names = F)
    sink()
    close(gz1)
}
q()
