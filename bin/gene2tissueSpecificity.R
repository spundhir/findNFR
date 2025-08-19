#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## parse command line arguments
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
option_list <- list(
	make_option(c("-g", "--genome"), default="mm10", help="genome (mm10 or hg38; default=%default)")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## check, if all required arguments are given
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
if(is.null(opt$genome)) {
	cat("\nProgram: gene2tissueSpecificity.R (R script to compute tissue specificity of a gene expression)\n")
	cat("Author: BRIC, University of Copenhagen, Denmark\n")
	cat("Version: 1.0\n")
	cat("Contact: pundhir@binf.ku.dk\n");
	print_help(parser)
	q()
}

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## load packages
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
suppressPackageStartupMessages(library(BgeeDB))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(biomaRt))

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## load custom functions
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
##-----------------------------------------##
# Calculate tau index (cell type specificity)
##-----------------------------------------##
#Function require a vector with expression of one gene in different tissues.
#If expression for one tissue is not known, gene specificity for this gene is NA
#Minimum 2 tissues
matrix2tau <- function(x)
{
  if(all(!is.na(x)))
  {
    if(min(x, na.rm=TRUE) >= 0)
    {
      if(max(x)!=0)
      {
        x <- (1-(x/max(x)))
        res <- sum(x, na.rm=TRUE)
        res <- res/(length(x)-1)
      } else {
        res <- 0
      }
    } else {
      res <- NA
      #print("Expression values have to be positive!")
    } 
  } else {
    res <- NA
    #print("No data for this gene avalable.")
  } 
  return(res)
}

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## start analysis based on input genome
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
CURRENT_DIR=getwd()
setwd("/home/xfd783/data/09_ALL_PUBLIC/database/bgee/rnaseq")
if(opt$genome=="mouse" | opt$genome=="mm10") {
  bgee <- Bgee$new(species = "Mus_musculus", dataType = "rna_seq", pathToData = "/home/xfd783/data/09_ALL_PUBLIC/database/bgee/rnaseq")
  data <- getData(bgee)
  expr <- dcast(as.data.table(data[,c("Gene.ID", "Anatomical.entity.name", "FPKM")]), Gene.ID ~ Anatomical.entity.name, value.var = "FPKM",
               fun.aggregate = function(x) mean(x, na.rm = TRUE))
  colnames(expr) <- gsub("[^A-Za-z0-9]+", "_", gsub("[^A-Za-z0-9]+$", "", gsub("\"", "", colnames(expr)))) %>% tolower
  mart <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl", version = 102)
  coor <- getBM(filters="ensembl_gene_id", 
                attributes=c("ensembl_gene_id", "external_gene_name", "chromosome_name", "start_position", "end_position", "strand", "gene_biotype"), 
                values=expr$gene_id, mart=mart) %>% mutate(chromosome_name = sprintf("chr%s", chromosome_name)) %>%
    'colnames<-'(gsub("_position|omosome_name", "", colnames(.))) %>% mutate(strand = ifelse(strand==1, "+", "-"))
  expr <- merge(coor, expr, by.x="ensembl_gene_id", by.y="gene_id", all.y=T)
  expr <- expr[which(expr$gene_biotype=="protein_coding"),]
  mat <- expr[,c(8:ncol(expr))]
  # mat <- expr[,c(grep(paste(c("cerebellum", "cerebral_cortex", "heart", "right_kidney", "liver", "lung", "placenta", "intestine", "spleen", "testis", "thymus", "adrenal", "bladder", "colon", "duodenum", "flobe", "gfat", "lgintestine", "mamgland", "ovary", "sfat", "stomach"), collapse = "|"), colnames(expr), value=T))]
  expr$tau <- apply(log2(mat+1), 1, matrix2tau)
  # hist(expr$tau)
} else if(opt$genome=="human" | opt$genome=="hg38") {
  bgee <- Bgee$new(species = "Homo_sapiens", dataType = "rna_seq", pathToData = "/home/xfd783/data/09_ALL_PUBLIC/database/bgee/rnaseq")
  data <- getData(bgee)
  expr <- dcast(as.data.table(data[,c("Gene.ID", "Anatomical.entity.name", "FPKM")]), Gene.ID ~ Anatomical.entity.name, value.var = "FPKM",
               fun.aggregate = function(x) mean(x, na.rm = TRUE))
  colnames(expr) <- gsub("[^A-Za-z0-9]+", "_", gsub("[^A-Za-z0-9]+$", "", gsub("\"", "", colnames(expr)))) %>% tolower
  mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = 114)
  coor <- getBM(filters="ensembl_gene_id", 
                attributes=c("ensembl_gene_id", "external_gene_name", "chromosome_name", "start_position", "end_position", "strand", "gene_biotype"), 
                values=expr$gene_id, mart=mart) %>% mutate(chromosome_name = sprintf("chr%s", chromosome_name)) %>%
    'colnames<-'(gsub("_position|omosome_name", "", colnames(.))) %>% mutate(strand = ifelse(strand==1, "+", "-"))
  expr <- merge(coor, expr, by.x="ensembl_gene_id", by.y="gene_id", all.y=T)
  expr <- expr[which(expr$gene_biotype=="protein_coding"),]
  mat <- expr[,c(8:ncol(expr))]
  # mat <- expr[,c(grep(paste(c("cerebellum", "cerebral_cortex", "heart", "right_kidney", "liver", "lung", "placenta", "intestine", "spleen", "testis", "thymus", "adrenal", "bladder", "colon", "duodenum", "flobe", "gfat", "lgintestine", "mamgland", "ovary", "sfat", "stomach"), collapse = "|"), colnames(expr), value=T))]
  expr$tau <- apply(log2(mat+1), 1, matrix2tau)
  # hist(expr$tau)
}
expr[expr == ""] <- NA
setwd(CURRENT_DIR)

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## output file containing gene expression specificity information
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
write.table(expr[,c(3:5,2,ncol(expr),6,7,1,8:(ncol(expr)-1))], file=sprintf("%s_tissueSpecificity.bed", opt$genome), sep="\t", col.names = T, row.names = F, quote = F)

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## save objects to .RDS file
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
saveRDS(list(expr=expr), file=sprintf("%s_tissueSpecificity.Rds", opt$genome))

q()
