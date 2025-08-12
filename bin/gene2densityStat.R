#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## parse command line arguments
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
option_list <- list(
	make_option(c("-i", "--inFile"), help="input file containing gene coordinates (can be stdin)"),
	make_option(c("-g", "--genome"), default="mm10", help="genome (mm10 or hg38; default=%default)")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## check, if all required arguments are given
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
if(is.null(opt$inFile)) {
	cat("\nProgram: gene2densityStat.R (R script to compute gene density in the genome)\n")
	cat("Author: BRIC, University of Copenhagen, Denmark\n")
	cat("Version: 1.0\n")
	cat("Contact: pundhir@binf.ku.dk\n");
	print_help(parser)
	q()
}

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## load packages
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(Hmisc))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(data.table))

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## step-1 initialize input files based on genome
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
if(identical(opt$inFile, "stdin")==T) {
    genes <- read.table(file("stdin"))[,c(1:6)] %>% 'colnames<-'(c("chr", "start", "end", "name", "score", "strand"))
} else {
    #GENE_FILE <- system(sprintf("source ~/.bashrc && initialize_genome -g %s -e", opt$genome), intern=T)
    #genes <- read.table(pipe(sprintf("grep -w protein_coding %s | grep -wv chrMT", GENE_FILE)))[,c(1:6)] %>% 'colnames<-'(c("chr", "start", "end", "name", "score", "strand"))
    genes <- read.table(opt$inFile)[,c(1:6)] %>% 'colnames<-'(c("chr", "start", "end", "name", "score", "strand"))
}
genes <- genes[match(unique(genes$name), genes$name),]
genes <- genes[order(genes$chr, genes$start, genes$end),]
GENOME_FILE <- system(sprintf("source ~/.bashrc && initialize_genome -g %s", opt$genome), intern=T)
TMP_FILE <- paste0("/tmp/file_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
write.table(genes, TMP_FILE, sep = "\t", quote = F, row.names = F, col.names = F)
## -k parameter (grep -w protein_coding /home/xfd783/genomes/annotations/BED/mm10_ensembl_gene.bed | cut -f 1 | sort | uniq -c | sed -E 's/^\s+//g' | sort -k 1rn,1)
## -k parameter can also be sort(table(genes$chr), decreasing = T)[1]
CMD <- sprintf("source ~/.bashrc && closestBed -a %s -b %s -d -k 2500 -N -g %s", TMP_FILE, TMP_FILE, GENOME_FILE)
df <- read.delim(text=system(CMD, intern = T), sep="\t", header = F)
df <- df[,c(1,4,10,13)] %>% 'colnames<-'(c("chr", "qGene", "sGene", "dist"))
invisible(suppressWarnings(file.remove(TMP_FILE)))

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## step-2 start analysis
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
p <- lapply(unique(df$chr), function(CHR) {
  df.tmp <- df[which(df$chr==CHR),c(2,3,4)]
  mat <- as.data.frame(dcast(as.data.table(df.tmp, 2500), qGene ~ sGene, value.var="dist"))
  row.names(mat) <- mat$qGene
  mat$qGene <- NULL
  mat[is.na(mat)] <- 0
  # head(apply((mat %>% 'diag<-'(NA)), 1, min, na.rm=T))
  cbind(apply((mat %>% 'diag<-'(NA)), 1, function(x) {
                minD <- min(x, na.rm=T)
              }) %>% as.data.frame,
        apply((mat %>% 'diag<-'(NA)), 1, function(x) {
                meanD = (mean(sort(unlist(x)) %>% head(10))) ## IMP PARAM
                minD <- min(x, na.rm=T)
                ((minD+1) * meanD) %>% log1p
              }) %>% as.data.frame
  ) %>% 'colnames<-'(c("dist2ClosestGene", "geneDensityScore"))
  # plot(apply((mat %>% 'diag<-'(NA)), 1, min, na.rm=T) %>% log1p,
  #      # apply(mat, 1, function(x) (mean(sort(unlist(x)) %>% head(10)+1) %>% log))
  #      log1p((apply((mat %>% 'diag<-'(NA)), 1, min, na.rm=T)+1)*apply(mat, 1, function(x) (mean(sort(unlist(x)) %>% head(10)))))
  #      )
})
df.res <- do.call(rbind, p) %>% mutate(gene=row.names(.)) %>% dplyr::select(c("gene", "dist2ClosestGene", "geneDensityScore"))
df.res <- merge(genes, df.res, by.x="name", by.y="gene")
df.res <- df.res[,c(2:4,1,5:ncol(df.res))]
df.res$geneLength <- df.res$end - df.res$start
df.res$geneDensityClass <- cut2(df.res$geneDensityScore, g=10, levels.mean=T)
levels(df.res$geneDensityClass) <- cut2(df.res$geneDensityScore, g=10, onlycuts = T)[2:length(cut2(df.res$geneDensityScore, g=10, onlycuts = T))]
levels(df.res$geneDensityClass) <- seq(length(levels(df.res$geneDensityClass)), 1)
# plot(log(df.res$dist2ClosestGene+1), df.res$geneDensityScore, xlab="Distance to closest gene (log)", ylab="Gene density score (norm)")
# ggdensity(df.res, x="dist2ClosestGene+1", xscale="log2", fill="#9ebcda", facet.by="chr")

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## step-3 output file containing gene density information
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
write.table(df.res, "", sep="\t", quote=F, row.names = F, col.names=T)

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## save objects to .RDS file
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
saveRDS(list(df=df, df.res=df.res), file=sprintf("%s_geneDensityStat.Rds", opt$genome))

q()
