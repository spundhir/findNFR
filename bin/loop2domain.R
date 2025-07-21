#!/usr/bin/env Rscript
#suppressPackageStartupMessages(library("optparse", lib.loc = "/home/users/sachin/R/x86_64-redhat-linux-gnu-library/2.15"))
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
  make_option(c("-i", "--inDomainFile"), help="input file containing gene interaction domains"),
  make_option(c("-g", "--genome"), default="mm10", help="genome for which to perform the analysis (mm10, hg38; default: mm10)"),
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
suppressPackageStartupMessages(library("AnnotationDbi"))

## read input file
if(identical(opt$inDistFile, "stdin")==T) {
    df <- read.table(file("stdin"), header=T)
} else {
    df <- read.table(opt$inDomainFile, header=T)
}
colnames(df) <- c("chr", "start", "end", "targetGene", "iScore", "strand", "sampleCount")

## determine gene order sorted by their genomic position
GENES <- (df %>% separate_longer_delim(targetGene, delim = ",") %>% dplyr::select(targetGene) %>% unique %>% unlist %>% unique)
if(opt$genome=="mm10") {
    #library(org.Mm.eg.db)
    #library(TxDb.Mmusculus.UCSC.mm10.knownGene)
    ## Convert gene symbols to Entrez IDs
    #entrez_ids <- mapIds(org.Mm.eg.db,
    #                        keys = GENES,
    #                        column = "ENTREZID",
    #                        keytype = "SYMBOL",
    #                        multiVals = "first") %>% as.data.frame %>% 'colnames<-'("entrez_id") %>% mutate(gene_name=row.names(.))
    #GENES_SORTED <- merge(genes(TxDb.Mmusculus.UCSC.mm10.knownGene, single.strand.genes.only=FALSE) %>% 
    #                            as.data.frame %>% 
    #                            dplyr::filter(group_name %in% entrez_ids$entrez_id), 
    #                        entrez_ids, by.x="group_name", by.y="entrez_id") %>% 
    #                    dplyr::arrange(seqnames, start, end) %>% dplyr::select(gene_name) %>% unlist
    GENES_SORTED <- read.table("/scratch/genomes/annotations/BED/mm10_ensembl_tss.bed") %>% dplyr::arrange(V1, V2, V3) %>% dplyr::select(V4) %>% unique %>% unlist
    GENES_SORTED <- GENES_SORTED[GENES_SORTED %in% GENES]
} else if(opt$genome=="hg38") {
    #library(org.Hs.eg.db)
    #library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    ## Convert gene symbols to Entrez IDs
    #entrez_ids <- mapIds(org.Hs.eg.db,
    #                        keys = GENES,
    #                        column = "ENTREZID",
    #                        keytype = "SYMBOL",
    #                        multiVals = "first") %>% as.data.frame %>% 'colnames<-'("entrez_id") %>% mutate(gene_name=row.names(.))
    #GENES_SORTED <- merge(genes(TxDb.Hsapiens.UCSC.hg38.knownGene, single.strand.genes.only=FALSE) %>% 
    #                            as.data.frame %>% 
    #                            dplyr::filter(group_name %in% entrez_ids$entrez_id), 
    #                        entrez_ids, by.x="group_name", by.y="entrez_id") %>% 
    #                    dplyr::arrange(seqnames, start, end) %>% dplyr::select(gene_name) %>% unlist
    GENES_SORTED <- read.table("/scratch/genomes/annotations/BED/hg38_ensembl_tss.bed") %>% dplyr::arrange(V1, V2, V3) %>% dplyr::select(V4) %>% unique %>% unlist
    GENES_SORTED <- GENES_SORTED[GENES_SORTED %in% GENES]
} else {
    cat("Unknown genome provided\n");
    print_help(parser)
    q()
}

## organize the gene co-interaction matrix
if(opt$processors==1) {
    dl <- lapply(GENES, function(GENE) {
        df[grepl(sprintf("(^|,)%s(,|$)", GENE), df$targetGene),c("targetGene","iScore")] %>% separate_longer_delim(cols=c(targetGene,iScore), delim = ",") %>% mutate(iScore=as.numeric(iScore)) %>% summaryBy(formula = iScore ~ targetGene, FUN = "sum") %>% mutate(gene=GENE) %>% 'colnames<-'(c("sGene", "iScore", "qGene")) %>% dplyr::select(c("qGene", "sGene", "iScore"))
      })
} else {
    dl <- mclapply(GENES, function(GENE) {
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
write.table(mat[GENES_SORTED,GENES_SORTED], file="", col.names=T, row.names=T, quote=F, sep="\t")

q()
