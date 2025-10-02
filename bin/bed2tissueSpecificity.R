#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## parse command line arguments
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
option_list <- list(
	make_option(c("-i", "--inFile"), help="bed file (can be stdin) [format: chr start end name score (comma sep.) strand tissue (comma sep.)]"),
	make_option(c("-M", "--useMeuleman"), action="store_true", help="compute using data from Meuleman et al. 2020")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## check, if all required arguments are given
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
if(is.null(opt$inFile)) {
	cat("\nProgram: bed2tissueSpecificity.R (R script to compute tissue specificity using score and tissue information)\n")
	cat("Author: BRIC, University of Copenhagen, Denmark\n")
	cat("Version: 1.0\n")
	cat("Contact: pundhir@binf.ku.dk\n");
	print_help(parser)
	q()
}

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## load packages
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyr))

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## load custom functions
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
##-----------------------------------------##
## Calculate tau index (cell type specificity)
## Paper: https://academic.oup.com/bib/article/18/2/205/2562739
##-----------------------------------------##
#Function require a vector with expression of one gene in different tissues.
#If expression for one tissue is not known, gene specificity for this gene is NA
#Minimum 2 tissues
vector2tau <- function(x)
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

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## read input file
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
if (identical(opt$inFile, "stdin") == T) {
  df <- read.table(file("stdin"), header = F)
} else {
  df <- read.table(opt$inFile, header = F)
}
## check if file has header (all elements in the first row are non-numeric)
if(length(grep("FALSE", unlist(lapply(df[1,], function(x) is.na(suppressWarnings(as.numeric(x)))))))==0) {
    colnames(df) <- unlist(df[1,]); 
    df <- df[-1,];
} else {
    df <- df[,c(1:7)]
    colnames(df) <- c("chr", "start", "end", "name", "score", "strand", "tissue")
}
#df <- read.table(pipe("zcat /home/xfd783/data/09_ALL_PUBLIC/database/02_accessibilityCatalog/chipAtlas/hg38/hg38_atac_tissueSpecificity.bed.gz | head -n 10"), header=T)[,c(1:7)] %>% 'colnames<-'(c("chr", "start", "end", "name", "score", "strand", "tissue"))

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## start analysis based on input bed file
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
if(is.null(opt$useMeuleman)) {
    dfE <- cbind((df %>% separate_longer_delim(score, delim = ",") %>% dplyr::select(!tissue)), (df %>% separate_longer_delim(tissue, delim = ",") %>% dplyr::select(tissue)))
    dfE <- aggregate(as.numeric(score) ~ chr + start + end + name + strand + tissue, dfE, FUN = mean) %>% 
              'colnames<-'(c("chr", "start", "end", "name", "strand", "tissue", "score")) %>% 
              dplyr::select(c("chr", "start", "end", "name", "score", "strand", "tissue")) %>% 
              dplyr::arrange(chr, start, end, name)
    mat <- matrix(0, length(unique(dfE$name)), length(unique(dfE$tissue)), dimnames = list(unique(dfE$name), unique(dfE$tissue)))
    mat[cbind(dfE$name, dfE$tissue)] <- as.numeric(dfE$score)
    mat <- as.data.frame(mat)
    mat$tau <- apply(log2(mat+1), 1, vector2tau)
    mat$meanScore <- rowMeans(mat)
    df <- merge(df, mat, by.x="name", by.y="row.names")
    df <- df[,c(2:4,1,5:ncol(df))]
} else {
    counts <- fread("~/data/09_ALL_PUBLIC/database/02_accessibilityCatalog/dhsCatalog/hg38/dat_FDR01_hg38.txt.gz", header=F)
    metaData <- fread("~/data/09_ALL_PUBLIC/database/02_accessibilityCatalog/dhsCatalog/hg38/DHS_Index_and_Vocabulary_metadata.tsv", header=T)
    
    mat <- as.data.frame(lapply(unique(metaData[which(metaData$Organ!=""),]$Organ), function(x) {
                                                                          i <- which(metaData$Organ==x)
                                                                          rowMeans(counts[,..i])
                                                                        }
                                )
                         )
    df$tau <- apply(mat, 1, vector2tau)
    df$gini <- apply(mat, 1, function(x) DescTools::Gini(x, na.rm = T))
    # ggscatter_with_contour(df, var_x="numsamples", var_y="tau")
}

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## output file containing tissue specificity information
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
write.table(df, file="", sep="\t", col.names = T, row.names = F, quote = F)

q()
