#!/usr/bin/env Rscript
#suppressPackageStartupMessages(library("optparse", lib.loc = "/home/users/sachin/R/x86_64-redhat-linux-gnu-library/2.15"))
suppressPackageStartupMessages(library("optparse"))

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## parse command line arguments
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
option_list <- list(
  make_option(c("-i", "--inFile"), help="input file containing gene interaction info (chr start end target_gene(s) iScore(s) strand #sample(s)"),
  make_option(c("-g", "--genome"), default="mm10", help="genome for which to perform the analysis (mm10, hg38; default: mm10)"),
  make_option(c("-p", "--processors"), default=1, help="number of processors to use (default: %default)"),
  make_option(c("-n", "--name"), help="string, which will be used to name output files (default: inFile)")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## check, if all required arguments are given
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
if(is.null(opt$inFile)) {
  cat("\nProgram: loop2domain.R (organize gene x gene interaction domain matrix)\n")
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
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("doBy"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("AnnotationDbi"))
suppressPackageStartupMessages(library("apcluster"))
suppressPackageStartupMessages(library("lsa"))
suppressPackageStartupMessages(library("GenomicRanges"))

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## load custom functions
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
bed2overlap <- function(df, bed_file, win_flank=0, minOverlap=1, selectFirstOverlap=F) {
  ## df: input data frame (format: chr start end unique_id) OR promoter_cor (chr:start-end) unique_id)
  ## bed_file: region to check for overlap; can also be a data frame (format: chr start end)
  if(ncol(df)==2) {
    t <- as.data.frame(t(as.data.frame(lapply(df[,1], function(x)  unlist(strsplit(x, "[:-]+"))))))
    row.names(t) <- row.names(df)
    t$name <-   df$name
    df <- t
    df$V2 <- as.numeric(df$V2)
    df$V3 <- as.numeric(df$V3)
  }
  
  # df <- df[,c(1:4)]
  colnames(df)[1:4] <- c("chr", "start", "end", "name")
  colnames(df) <- sprintf("%s_q", colnames(df))
  if(win_flank>0) {
    df$start_q <- (df$start_q+df$end_q)/2
    df$end_q <- df$start_q+win_flank
    df$start_q <- df$start_q-win_flank
  }
  
  if(!is.data.frame(bed_file)) {
    bed_file <- read.table(pipe(sprintf("grep -vi start %s", bed_file)), header=F, sep="\t")
  }
  
  # bed_file <- bed_file[,c(1:3)]
  colnames(bed_file)[1:3] <- c("chr", "start", "end")
  colnames(bed_file) <- sprintf("%s_s", colnames(bed_file))
  
  ## which(df$chr_q %in% unique(df[which(df$chr_q %in% bed_file$chr_s),]$chr_q)) ## no need of this
  query = makeGRangesFromDataFrame(df[,c("chr_q", "start_q", "end_q", "name_q")],
                                   seqnames.field = "chr_q", start.field = "start_q", end.field = "end_q",
                                   keep.extra.columns = T, ignore.strand = T)
  
  subject = makeGRangesFromDataFrame(bed_file,
                                     seqnames.field = "chr_s", start.field = "start_s", end.field = "end_s",
                                     keep.extra.columns = T)
  
  if(selectFirstOverlap==T) {
    t <- data.frame(findOverlaps(query, subject, ignore.strand=T, select="first", minoverlap = minOverlap))
    t <- t %>% mutate(queryHits=row.names(t), subjectHits=t[,1])
    t <- t[which(!is.na(t[,1])),c("queryHits", "subjectHits")]
  } else {
    t <- data.frame(findOverlaps(query, subject, ignore.strand=T, select="all", minoverlap = minOverlap))
  }
  
  return(cbind(df[t$queryHits,], bed_file[t$subjectHits,]))
}

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## read input file
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
if(identical(opt$inDistFile, "stdin")==T) {
    df <- read.table(file("stdin"), header=T)
} else {
    df <- read.table(opt$inFile, header=T)
}
colnames(df) <- c("chr", "start", "end", "targetGene", "iScore", "strand", "sampleCount")

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## define output file name
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
if(is.null(opt$name)) {
    opt$name <- opt$inFile
}

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## determine gene order sorted by their genomic position
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
#GENES <- (df %>% separate_longer_delim(targetGene, delim = ",") %>% dplyr::select(targetGene) %>% unique %>% unlist %>% unique)
GENES <- unique(df$targetGene)
if(opt$genome=="mm10") {
    #library(org.Mm.eg.db)
    #library(TxDb.Mmusculus.UCSC.mm10.knownGene)
    ## Convert gene symbols to Entrez IDs
    #entrez_ids <- mapIds(org.Mm.eg.db,
    #                        keys = GENES,
    #                        column = "ENTREZID",
    #                        keytype = "SYMBOL",
    #                        multiVals = "first") %>% as.data.frame %>% 'colnames<-'("entrez_id") %>% mutate(gene_name=row.names(.))
    #GENE_SORTED <- merge(genes(TxDb.Mmusculus.UCSC.mm10.knownGene, single.strand.genes.only=FALSE) %>% 
    #                            as.data.frame %>% 
    #                            dplyr::filter(group_name %in% entrez_ids$entrez_id), 
    #                        entrez_ids, by.x="group_name", by.y="entrez_id") %>% 
    #                    dplyr::arrange(seqnames, start, end) %>% dplyr::select(gene_name) %>% unlist
    GENE_COOR <- read.table("/scratch/genomes/annotations/BED/mm10_ensembl_gene.bed")[,c(1:6)] %>% 'colnames<-'(c("chr", "start", "end", "name", "score", "strand")) %>% arrange(name) %>% distinct(name, .keep_all=T)
    GENE_SORTED <- read.table("/scratch/genomes/annotations/BED/mm10_ensembl_tss.bed") %>% dplyr::arrange(V1, V2, V3) %>% dplyr::select(V4) %>% unique %>% unlist
    GENE_SORTED <- GENE_SORTED[GENE_SORTED %in% GENES]
} else if(opt$genome=="hg38") {
    #library(org.Hs.eg.db)
    #library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    ## Convert gene symbols to Entrez IDs
    #entrez_ids <- mapIds(org.Hs.eg.db,
    #                        keys = GENES,
    #                        column = "ENTREZID",
    #                        keytype = "SYMBOL",
    #                        multiVals = "first") %>% as.data.frame %>% 'colnames<-'("entrez_id") %>% mutate(gene_name=row.names(.))
    #GENE_SORTED <- merge(genes(TxDb.Hsapiens.UCSC.hg38.knownGene, single.strand.genes.only=FALSE) %>% 
    #                            as.data.frame %>% 
    #                            dplyr::filter(group_name %in% entrez_ids$entrez_id), 
    #                        entrez_ids, by.x="group_name", by.y="entrez_id") %>% 
    #                    dplyr::arrange(seqnames, start, end) %>% dplyr::select(gene_name) %>% unlist
    GENE_COOR <- read.table("/scratch/genomes/annotations/BED/hg38_ensembl_gene.bed")[,c(1:6)] %>% 'colnames<-'(c("chr", "start", "end", "name", "score", "strand")) %>% arrange(name) %>% distinct(name, .keep_all=T)
    GENE_SORTED <- read.table("/scratch/genomes/annotations/BED/hg38_ensembl_tss.bed") %>% dplyr::arrange(V1, V2, V3) %>% dplyr::select(V4) %>% unique %>% unlist
    GENE_SORTED <- GENE_SORTED[GENE_SORTED %in% GENES]
} else {
    cat("Unknown genome provided\n");
    print_help(parser)
    q()
}

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## organize the gene co-interaction matrix
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
dfO <- bed2overlap(df, df, minOverlap = 5000)

dl <- mclapply(GENES, function(GENE) {
    #df[grepl(sprintf("(^|,)%s(,|$)", GENE), df$targetGene),c("targetGene","iScore")] %>% separate_longer_delim(cols=c(targetGene,iScore), delim = ",") %>% mutate(iScore=as.numeric(iScore)) %>% summaryBy(formula = iScore ~ targetGene, FUN = "sum") %>% mutate(gene=GENE) %>% 'colnames<-'(c("sGene", "iScore", "qGene")) %>% dplyr::select(c("qGene", "sGene", "iScore"))
  dfO[which(dfO$name_q==GENE),] %>% dplyr::mutate(iScore=log(iScore_q * iScore_s)) %>% dplyr::select(targetGene_s, iScore) %>% summaryBy(formula = iScore ~ targetGene_s, FUN = "sum") %>% mutate(gene=GENE) %>% 'colnames<-'(c("sGene", "iScore", "qGene")) %>% dplyr::select(c("qGene", "sGene", "iScore"))
  }, mc.cores=opt$processors, mc.set.seed = 5)

## Combine all into one table
dt <- rbindlist(dl)

## Get unique genes
all_genes <- unique(dt$qGene)

# Initialize matrix
mat <- matrix(NA, length(all_genes), length(all_genes), dimnames = list(all_genes, all_genes))

# Fill
mat[cbind(dt$qGene, dt$sGene)] <- dt$iScore

mat <- mat[GENE_SORTED,GENE_SORTED]
mat.norm <- sweep(mat, 1, diag(as.matrix(mat)), "/")

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## write matrix to file
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
write.table(mat, file=sprintf("%s.matrix", opt$name), col.names=T, row.names=T, quote=F, sep="\t")

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## determine domains
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
# gene_chr <- (df %>% separate_longer_delim(targetGene, delim = ",") %>% dplyr::select(chr, targetGene) %>% unique)
gene_chr <- (df %>% dplyr::select(chr, targetGene) %>% unique)

t <- mclapply(unique(gene_chr$chr), function(CHR) {
  mat.sub <- mat.norm[which(colnames(mat.norm) %in% gene_chr[which(gene_chr==CHR),]$targetGene), 
                 which(colnames(mat.norm) %in% gene_chr[which(gene_chr==CHR),]$targetGene)] %>%
    replace(is.na(.), 0)
  
  ## APPROACH-1
  # t <- lsa::cosine(t(mat.sub))
  # t <- lsa::cosine(t(matrix2norm(t(mat.sub), method = "softmax")))
  # pheatmap(t, cluster_rows = F, cluster_cols = F, border_color = NA, scale="none")
  # CI <- cluster_matrix(t, y = 114, usePheatmap = T)
  # CI$cluster[which(names(CI$cluster)=="Meis1")]
  # which(CI$cluster==CI$cluster[which(names(CI$cluster)=="Meis1")])
  # barplot(CI$withinss, names=unique(sort(CI$cluster)), las=2)
  # barplot(CI$centers[1,])
  # ggarrange(matrix2Heatmap(mat.sub, y = CI$cluster, clusterRows = F, clusterCols = F))
  
  ## APPROACH-2
  # mat.sub <- mat.sub[row.names(mat.sub[rowSums(mat.sub != 0) > 1, ]),row.names(mat.sub[rowSums(mat.sub != 0) > 1, ])]
  # ap_result <- apcluster(s = as.matrix(t(mat.sub)), seed=10, maxits=10000)
  ap_result <- apcluster(cosine, t(mat.sub), seed=10, maxits=10000, details=TRUE)
  geneCluster <- do.call(rbind, lapply(seq(1:length(ap_result@clusters)), function(i) cbind(names(ap_result@clusters[[i]]), i))) %>% as.data.frame %>% 'colnames<-'(c("gene", "cluster"))
  geneBlackListed <- lapply(seq(1:length(ap_result@clusters)), function(i) unlist(lapply(names(ap_result@clusters[[i]]), function(x) { if(length(which(mat.sub[x, names(ap_result@clusters[[i]])]>0))==1) { x; } }))) %>% unlist
  if(length(geneBlackListed)>0) {
    geneCluster[which(geneCluster$gene %in% geneBlackListed),]$cluster <- 0
  }
  geneCluster <- geneCluster[order(geneCluster$cluster),]
  geneCluster$cluster <- sprintf("%s_%s", geneCluster$cluster, CHR)
  geneCluster
  ## sanity check
  # GENE="Hoxa9"
  # geneCluster[which(geneCluster$gene==GENE),]$cluster
  # geneCluster[which(geneCluster$cluster==geneCluster[which(geneCluster$gene==GENE),]$cluster),]$gene
  # mat.sub[GENE, geneCluster[which(geneCluster$cluster==geneCluster[which(geneCluster$gene==GENE),]$cluster),]$gene]
  # mat.sub[GENE,][which(mat.sub[GENE,]>0)]
  # ggarrange(matrix2Heatmap(mat.sub, y=geneCluster$cluster, clusterRows = F, clusterCols = F))
  # pheatmap(as.matrix(mat.sub[geneBlackListed,]), cluster_rows = F, cluster_cols = F, border_color = NA)
  
  ## APPROACH-3
  # umap_result <- uwot::umap(t, n_neighbors = 15, min_dist = 0.1, metric = "euclidean")
  # # Put in data frame for plotting
  # umap_df <- data.frame(UMAP1 = umap_result[,1],
  #                       UMAP2 = umap_result[,2],
  #                       Gene = rownames(mat.sub))
  # ggscatter(umap_df, x="UMAP1", y="UMAP2", label = "Gene", label.select = "Meis1")
}, mc.cores=opt$processors, mc.set.seed = 5)
geneCluster <- do.call(rbind, t)
# (geneCluster %>% dplyr::filter(cluster==(geneCluster %>% dplyr::filter(gene=="Myc") %>% dplyr::select(cluster))$cluster))$gene
geneCluster <- merge(GENE_COOR, geneCluster, by.x="name", by.y="gene") %>% dplyr::select(c("chr", "start", "end", "name", "score", "strand", "cluster"))
geneCluster <- merge(geneCluster, diag(as.matrix(mat)) %>% as.data.frame %>% 'colnames<-'(c("iScore")), by.x="name", by.y="row.names") %>% dplyr::select(c("chr", "start", "end", "name", "iScore", "strand", "cluster"))

geneCluster <- geneCluster %>% group_by(cluster) %>% summarise(chr = unique(chr, na.rm = TRUE),
                                                start = min(start, na.rm = TRUE),
                                                end = max(end, na.rm = TRUE),
                                                gene = paste(name, collapse = ","),
                                                iScore = paste(sprintf("%0.2f", iScore), collapse = ","),
                                                strand = "."
                                                ) %>% dplyr::select(c("chr", "start", "end", "gene", "iScore", "strand", "cluster")) %>% dplyr::arrange(chr, start, end)

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## write gene domains to file
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
write.table(geneCluster, file=sprintf("%s.domain", opt$name), col.names=T, row.names=F, quote=F, sep="\t")

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## save objects to .RDS file
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
saveRDS(list(df=df, mat=mat, gene_chr=gene_chr, geneCluster=geneCluster), file=sprintf("%s.Rds", opt$name))

q()
