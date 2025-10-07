#!/usr/bin/env Rscript
#suppressPackageStartupMessages(library("optparse", lib.loc = "/home/users/sachin/R/x86_64-redhat-linux-gnu-library/2.15"))
suppressPackageStartupMessages(library("optparse"))

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## parse command line arguments
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
option_list <- list(
  make_option(c("-i", "--inFile"), help="input file containing peaks in bed format [format: chr start end name score strand]"),
  make_option(c("-g", "--genome"), default="mm10", help="genome (mm10 or hg38; default=%default)")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## check, if all required arguments are given
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
if(is.null(opt$inFile)) {
  cat("\nProgram: peaks2spatialAnnotation.R (annotate peaks classified based on distance to TSS)\n")
  cat("Author: BRIC, University of Copenhagen, Denmark\n")
  cat("Version: 1.0\n")
  cat("Contact: pundhir@binf.ku.dk\n");
  print_help(parser)
  q()
}

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## load packages
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("ggrastr"))
suppressPackageStartupMessages(library("ggpubr"))
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("Hmisc"))
suppressPackageStartupMessages(library("cowplot"))

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## load custom functions
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
suppressPackageStartupMessages(source("~/software/myScripts/rScripts/FUNCTIONS.R"))

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## initialize input files based on genome
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## read input file
if(identical(opt$inFile, "stdin")==T) {
    peaks <- read.table(file("stdin"))
} else {
    peaks <- read.table(opt$inFile)
}
# peaks <- read.table(pipe("zcat ~/data/09_ALL_PUBLIC/database/02_accessibilityCatalog/dhsCatalog/mm10_ocr.bed.gz | head -n 1000"), header=F)

## check if file has header (all elements in the first row are non-numeric)
if(length(grep("FALSE", unlist(lapply(peaks[1,], function(x) is.na(suppressWarnings(as.numeric(x)))))))==0) {
  colnames(peaks) <- peaks[1,]
  peaks <- peaks[-1,];
} else {
  colnames(peaks) <- unlist(lapply(seq(1,ncol(peaks)), function(x) sprintf("description_%d",x)))
}
colnames(peaks)[c(1:6)] <- c("chr", "start", "end", "name", "score", "strand")

GENOME_FILE <- system(sprintf("source ~/.bashrc && initialize_genome -g %s", opt$genome), intern=T)
TSS_FILE <- system(sprintf("source ~/.bashrc && initialize_genome -g %s -t", opt$genome), intern=T)
GENE_TAU_FILE <- system(sprintf("source ~/.bashrc && initialize_genome -g %s -S", opt$genome), intern=T)
DHS_TAU_FILE <- system(sprintf("source ~/.bashrc && initialize_genome -g %s -H", opt$genome), intern=T)
GENEDENSITY_FILE <- system(sprintf("source ~/.bashrc && initialize_genome -g %s -M", opt$genome), intern=T)
CPG_FILE <- system(sprintf("source ~/.bashrc && initialize_genome -g %s -c", opt$genome), intern=T)

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## start analysis
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## closest/target gene information
df <- merge(peaks, linkDHS2Genes(bed2window(peaks, win = 250, flank_to_tss = F), genome = opt$genome, useLoops = F, strandAware = T)[,c("name", "closest_gene", "dist_to_closest_gene")],
            by.x="name", by.y="name")
df <- df[,c(2:4,1,5:ncol(df))]
df <- merge(df, linkDHS2Genes(bed2window(peaks, win = 250, flank_to_tss = F), genome = opt$genome, useLoops = T, minoverlap = 100)[,c("name", "target_gene", "cInteraction_score", "dist_to_target_gene")],
  by.x="name", by.y="name")
df <- df[,c(2:4,1,5:ncol(df))]
df$dist_to_closest_gene <- log(abs(df$dist_to_closest_gene)+1)
df <- merge((df %>% dplyr::select(!dist_to_target_gene)), (df %>% separate_longer_delim(dist_to_target_gene, ",") %>% mutate(dist_to_target_gene = log(as.numeric(dist_to_target_gene)+1)) %>%
                                                       aggregate(dist_to_target_gene ~ name, toString) %>% mutate(dist_to_target_gene = gsub("\\s+", "", dist_to_target_gene))),
            by.x="name", by.y="name", all.x=T)

## closest/target gene tissue specificity information
tmp <- read.table(GENE_TAU_FILE, header=T)[,c("external_gene_name", "tau")] %>% dplyr::filter(!is.na(external_gene_name))
df$closest_gene_tau <- tmp$tau[match(df$closest_gene, tmp$external_gene_name)]
t <- df %>% separate_longer_delim(target_gene, ",")
t$target_gene_tau <- round(tmp$tau[match(t$target_gene, tmp$external_gene_name)],7)
df <- merge(df, aggregate(target_gene_tau ~ name, t, toString) %>% mutate(target_gene_tau = gsub("\\s+", "", target_gene_tau)), by.x="name", by.y="name", all.x=T)
df <- df[,c(2:4,1,5:ncol(df))]

## DHS tissue specificity information
df <- merge(df, bed2overlap(df, bed_file = read.table(pipe(sprintf("zcat %s | cut -f 1-6", DHS_TAU_FILE)), header=T), selectFirstOverlap = T)[,c("name_q", "tau_s")],
            by.x="name", by.y="name_q", all.x=T) %>% rename(dhs_tau=tau_s)
df <- df[,c(2:4,1,5:ncol(df))]

## closest gene density information
df <- merge(df, read.table(GENEDENSITY_FILE, header=T)[,c("name", "geneDensityScore", "geneLength", "geneDensityClass")], by.x="closest_gene", by.y="name")
df$geneDensityClass <- factor(df$geneDensityClass)
df <- df[,c(2:8,1,9:ncol(df))]

## CpG island information
df$cpgOverlap <- ifelse(df$name %in% bed2overlap(df, CPG_FILE)$name_q, "CpG", "nonCpG")

## peak annoatation
df$annot.type <- "promoter"
df[which(abs(df$dist_to_closest_gene) > log(1000) & abs(df$dist_to_closest_gene) <= log(50000)),]$annot.type <- "proximal"
df[which(abs(df$dist_to_closest_gene) > log(50000)),]$annot.type <- "distal"
df$annot.type <- factor(df$annot.type, levels=c("promoter", "proximal", "distal"), ordered=T)
# table(df$annot.type)

## class defined based on distance to closest gene
# df$distClass <- cut2(abs(df$dist_to_closest_gene), cuts=log(c(0, 1000, 2500, 5000, 10000, 50000, 100000, 500000, 1000000, 40000000)), levels.mean = T)
df$distClass <- cut2(abs(df$dist_to_closest_gene), g=10, levels.mean = T)
levels(df$distClass) <- sprintf("%s_%s", seq(1,length(levels(df$distClass))), round(as.numeric(gsub("\\s+", "", levels(df$distClass))),2))
# table(df$distClass)

# ggscatter_with_contour(df, var_x="dhs_tau", var_y="closest_gene_tau", color_class = "distClass")
# ggboxplot(df, x="geneDensityClass", y="dhs_tau")
# ggboxplot(df %>% separate_longer_delim(target_gene_tau, ",") %>% mutate(target_gene_tau = as.numeric(target_gene_tau)), x="distClass", y="target_gene_tau")
# df %>% separate_longer_delim(dist_to_target_gene, ",") %>% 
#   mutate(dist_to_target_gene = log(as.numeric(dist_to_target_gene)+1)) %>% 
#   mutate(dist_to_closest_gene = (abs(dist_to_closest_gene))) %>% 
#   ggscatter(x="dist_to_target_gene", y="dist_to_closest_gene")

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## save output files
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
# ggsave(filename = opt$outFile, plot = P, width = 9, height = 6)
write.table(df, "", sep="\t", col.names = T, row.names = F, quote = F)
q()
