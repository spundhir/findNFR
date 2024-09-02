#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
	make_option(c("-i", "--inFile"), help="input file (format: target_gene [atac_signal] score class{1,}) (can be stdin)"),
	make_option(c("-C", "--classInfoCharacter"), help="class information is in character format", action="store_true"),
	make_option(c("-o", "--outRdsFile"), help="output .rds file to store results"),
	make_option(c("-g", "--organism"), help="organism for which to perform the analysis (mouse or human)"),
    make_option(c("-d", "--geneIdType"), default="SYMBOL", help="input gene id type (default=%default)"),
    make_option(c("-a", "--annotation"), default="C5", help="annotation type (default=%default)"),
    make_option(c("-p", "--pValue"), default="1e-05", help="p-value (defaut=%default)"),
	make_option(c("-c", "--processors"), default=1, help="number of processors to use"),
    make_option(c("-m", "--maxClass"), default=20, help="maximum number of go classes to plot (default=%default)"),
    make_option(c("-n", "--minGene"), default=10, help="minimum genes in a go class (default=%default)"),
    make_option(c("-w", "--figWidth"), default="10", help="width of output figure"),
    make_option(c("-t", "--figHeight"), default="20", help="height of output figure"),
	make_option(c("-r", "--ftrResFile"), help="input file containing manually filtered GO categories to plot"),
	make_option(c("-q", "--logCount"), help="plot count data in log", action="store_true")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if((is.null(opt$inFile) | is.null(opt$organism) | is.null(opt$outRdsFile))) {
	cat("\nProgram: multiClassGoAna.R (R script to compute GO analysis using GSEA)\n")
	cat("Author: BRIC, University of Copenhagen, Denmark\n")
	cat("Version: 1.0\n")
	cat("Contact: pundhir@binf.ku.dk\n");
	print_help(parser)
	q()
}

## load libraries
suppressPackageStartupMessages(library(DOSE))
suppressPackageStartupMessages(library(ReactomePA))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(mygene))
#suppressPackageStartupMessages(library(RDAVIDWebService)) #R CMD javareconf -e
suppressPackageStartupMessages(library(session))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(msigdbr))
suppressPackageStartupMessages(library(dplyr))
source("~/software/myScripts/R-scripts/FUNCTIONS.R")

## read input file
if(identical(opt$inFile, "stdin")==T) {
    df <- read.table(file("stdin"))
} else {
    df <- read.table(opt$inFile)
}

## format input data
if(opt$classInfoCharacter) {
    if(ncol(df)==3 & is.numeric(df[,2]) & is.numeric(df[,3])) {
        colnames(df) <- c("target_gene", "score", "class")
        t <- lapply(unique(df$class), function(x) {
            tmp <- df[which(!is.na(df$target_gene) & df$class==x),c("target_gene", "score", "class")] %>% separate_longer_delim(c(target_gene, score),  delim=",") %>% mutate(score=as.numeric(score)) %>% as.data.frame()
            tmp <- aggregate(score ~ target_gene, tmp, sum)
            colnames(tmp)[2] <- x
            return(tmp)
        })
        names(t) <- unique(df$class)
    } else if(ncol(df)==4 & is.numeric(df[,2]) & is.numeric(df[,3])) {
        colnames(df) <- c("target_gene", "signal", "score", "class")
        t <- lapply(unique(df$class), function(x) {
            tmp <- df[which(!is.na(df$target_gene) & df$class==x),c("target_gene", "signal", "score", "class")] %>% separate_longer_delim(c(target_gene, score),  delim=",") %>% mutate(score=as.numeric(score) * signal) %>% as.data.frame()
            tmp <- aggregate(score ~ target_gene, tmp, sum)
            colnames(tmp)[2] <- x
        return(tmp)
        })
        names(t) <- unique(df$class)
    } else if(ncol(df)>3 & is.numeric(df[,2]) & is.numeric(df[,3])) {
        colnames(df)[1:3] <- c("target_gene", "signal", "score")
        t <- lapply(colnames(df[4:ncol(df)]), function(x) {
            tmp <- df_enhancer[which(!is.na(df_enhancer$target_gene) & df_enhancer[,x]=="yes"),c("target_gene", "signal", "score")] %>% separate_longer_delim(c(target_gene, score),  delim=",") %>% mutate(class=x) %>% mutate(score=as.numeric(score) * signal) %>% as.data.frame()
            tmp <- aggregate(score ~ target_gene, tmp, sum)
            colnames(tmp)[2] <- x
            return(tmp)
        })
        names(t) <- colnames(df[4:ncol(df)])
    } else {
        colnames(df)[1:2] <- c("target_gene", "score")
        t <- lapply(colnames(df[3:ncol(df)]), function(x) {
            tmp <- df_enhancer[which(!is.na(df_enhancer$target_gene) & df_enhancer[,x]=="yes"),c("target_gene", "score")] %>% separate_longer_delim(c(target_gene, score),  delim=",") %>% mutate(class=x) %>% mutate(score=as.numeric(score)) %>% as.data.frame()
            tmp <- aggregate(score ~ target_gene, tmp, sum)
            colnames(tmp)[2] <- x
            return(tmp)
        })
        names(t) <- colnames(df[3:ncol(df)])
    }
    mat <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "target_gene", all.x = TRUE), t)[,-1]
    mat[is.na(mat)] <- 0
    mat <- as.data.frame(normalize.quantiles(as.matrix(mat)))
    row.names(mat) <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "target_gene", all.x = TRUE), t)$target_gene
    colnames(mat) <- colnames(Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "target_gene", all.x = TRUE), t)[,-1])
    mat <- log(mat+1)
    mat <- as.data.frame(scale(mat, center=T, scale=T))
} else {
    mat <- df[,c(-1)]
    row.names(mat) <- df[,1]
    mat[is.na(mat)] <- 0
}

## perform analysis
goRes <- multiClassGoAna(mat, org = opt$organism, caty = opt$annotation, sigT = opt$pValue, processors = opt$processors)
saveRDS(goRes, opt$outRdsFile)
q()
