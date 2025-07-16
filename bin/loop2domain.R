#!/usr/bin/env Rscript
#suppressPackageStartupMessages(library("optparse", lib.loc = "/home/users/sachin/R/x86_64-redhat-linux-gnu-library/2.15"))
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
  make_option(c("-i", "--inDomainFile"), help="input file containing gene interaction domains"),
  make_option(c("-p", "--processors"), default=1, help="number of processors to use (default: %default)")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if(is.null(opt$inDomainFile)) {
  cat("\nProgram: loop2domain.R (organize gene x gene interaction domain matrix)\n")
  cat("Author: BRIC, University of Copenhagen, Denmark\n")
  cat("Version: 1.0\n")
  cat("Contact: pundhir@binf.ku.dk\n");
  print_help(parser)
  q()
}

suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("ggrastr"))
suppressPackageStartupMessages(library("ggpubr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("doBy"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("parallel"))

## read input file
if(identical(opt$inDistFile, "stdin")==T) {
    df <- read.table(file("stdin"), header=T)
} else {
    df <- read.table(opt$inDomainFile, header=T)
}
colnames(df) <- c("chr", "start", "end", "targetGene", "iScore", "strand", "sampleCount")

## organize the gene co-interaction matrix
if(opt$processors==1) {
    dl <- lapply((df %>% separate_longer_delim(targetGene, delim = ",") %>% dplyr::select(targetGene) %>% unique %>% unlist %>% unique), function(GENE) {
        df[grepl(sprintf("(^|,)%s(,|$)", GENE), df$targetGene),c("targetGene","iScore")] %>% separate_longer_delim(cols=c(targetGene,iScore), delim = ",") %>% mutate(iScore=as.numeric(iScore)) %>% summaryBy(formula = iScore ~ targetGene, FUN = "sum") %>% mutate(gene=GENE) %>% 'colnames<-'(c("sGene", "iScore", "qGene")) %>% dplyr::select(c("qGene", "sGene", "iScore"))
      })
} else {
    dl <- mclapply((df %>% separate_longer_delim(targetGene, delim = ",") %>% dplyr::select(targetGene) %>% unique %>% unlist %>% unique), function(GENE) {
        df[grepl(sprintf("(^|,)%s(,|$)", GENE), df$targetGene),c("targetGene","iScore")] %>% separate_longer_delim(cols=c(targetGene,iScore), delim = ",") %>% mutate(iScore=as.numeric(iScore)) %>% summaryBy(formula = iScore ~ targetGene, FUN = "sum") %>% mutate(gene=GENE) %>% 'colnames<-'(c("sGene", "iScore", "qGene")) %>% dplyr::select(c("qGene", "sGene", "iScore"))
      }, mc.cores=opt$processors, mc.set.seed = 5)
}

## Combine all into one table
dt <- rbindlist(dl)

## Get unique genes
all_genes <- unique(dt$qGene)

# Initialize matrix
mat <- matrix(NA, length(all_genes), length(all_genes), dimnames = list(all_genes, all_genes))

# Fill
mat[cbind(dt$qGene, dt$sGene)] <- dt$iScore

## write matrix to stdout
write.table(mat, file="", col.names=T, row.names=T, quote=F, sep="\t")

q()
