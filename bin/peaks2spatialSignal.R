#!/usr/bin/env Rscript
#suppressPackageStartupMessages(library("optparse", lib.loc = "/home/users/sachin/R/x86_64-redhat-linux-gnu-library/2.15"))
suppressPackageStartupMessages(library("optparse"))

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## parse command line arguments
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
option_list <- list(
  make_option(c("-i", "--inFile"), help="input file containing peak score (macs2 output) [format: chr start end name score strand signalValue]"),
  make_option(c("-o", "--outFile"), help="output pdf file"),
  make_option(c("-g", "--genome"), default="mm10", help="genome (mm10 or hg38; default=%default)")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## check, if all required arguments are given
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
if(is.null(opt$inFile) | is.null(opt$outFile)) {
  cat("\nProgram: peaks2spatialSignal.R (plot peak signal classified based on distance to TSS)\n")
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
    peaks <- read.table(file("stdin"))[,c(1:7)]
} else {
    peaks <- read.table(opt$inFile)[,c(1:7)]
}
# peaks <- read.table("~/data/00_ALL_CHIP-SEQ_RAW/MLL-AF9/six1_on_peaks.bed", header=F)[,c(1:7)]

## check if file has header (all elements in the first row are non-numeric)
if(length(grep("FALSE", unlist(lapply(peaks[1,], function(x) is.na(suppressWarnings(as.numeric(x)))))))==0) {
  peaks <- peaks[-1,];
}
colnames(peaks) <- c("chr", "start", "end", "name", "score", "strand", "signalValue")

## reformat peak file to ensure correct format of signal values
if(length(which(!is.na(peaks$signalValue))) > 0 & length(grep("TRUE", is.na(suppressWarnings(as.numeric(peaks$signalValue)))))==0) {
  peaks$signalValue <- log2(as.numeric(peaks$signalValue)+1)
} else if(length(which(!is.na(peaks$score))) > 0 & length(grep("TRUE", is.na(suppressWarnings(as.numeric(peaks$score)))))==0) {
  peaks$signalValue <- log2(as.numeric(peaks$score)+1)
} else {
  cat("ERROR: no score column found\n");
  print_help(parser)
  q()
}

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
# df <- bed2closestCoor(peaks, read.table(pipe(sprintf("grep -w protein_coding %s | cut -f 1-6", TSS_FILE))), strandAware = T)[,c(1:7,11,14)]
# colnames(df) <- c("chr", "start", "end", "name", "score", "strand", "signalValue", "closest_gene", "dist_to_closest_gene")
df <- merge(peaks, linkDHS2Genes(bed2window(peaks, win = 250, flank_to_tss = F), genome = opt$genome, useLoops = F, strandAware = T)[,c("name", "closest_gene", "dist_to_closest_gene")],
            by.x="name", by.y="name")
df <- df[,c(2:4,1,5:ncol(df))]

df$dist_to_closest_gene <- log(abs(df$dist_to_closest_gene)+1) * ifelse(df$dist_to_closest_gene<0, -1, 1)
df$annot.type <- "promoter"
df[which(abs(df$dist_to_closest_gene) > log(1000) & abs(df$dist_to_closest_gene) <= log(50000)),]$annot.type <- "proximal"
df[which(abs(df$dist_to_closest_gene) > log(50000)),]$annot.type <- "distal"
df$annot.type <- factor(df$annot.type, levels=c("promoter", "proximal", "distal"), ordered=T)
# table(df$annot.type)

brk <- c(min(df$dist_to_closest_gene), -log(50000), -log(1000), 0, log(1000), log(50000), max(df$dist_to_closest_gene))
lbl_abs <- (table(cut(df$dist_to_closest_gene, 
                      breaks=c(min(df$dist_to_closest_gene), -log(50000), -log(1000), 0, log(1000), log(50000), max(df$dist_to_closest_gene)), 
                      labels=c(min(df$dist_to_closest_gene), -log(50000), 0, 0, log(50000), max(df$dist_to_closest_gene))))) %>% as.vector
lbl_per <- gsub("\\s", "", paste(sprintf("%0.1f", (lbl_abs*100)/nrow(df)), "%"))

## make the plot (location vs signal value)
MAX_VAL=round(max(df$signalValue),0)
p1 <- ggplot(df, aes(dist_to_closest_gene, signalValue)) +
            geom_point_rast(aes(color=annot.type), alpha=0.2) +
            stat_density_2d(geom = "polygon", contour = TRUE, aes(fill = after_stat(level)), colour = "black", alpha=0.5, show.legend = F) +
            scale_fill_distiller(palette = "Reds", direction = 1) + 
            #xlim(c(-15,15)) + 
            geom_vline(xintercept = brk, lty=2) +
            scale_x_continuous(breaks=brk, labels = as.numeric(sprintf("%0.0f", brk))) +
            #annotate("text", x = c(-10, -5, 0, 5, 10), y = 2, label = lbl,  parse = F, size=3) +  theme_classic2() +
            annotate("text", x = c(-12.5, -9, 0, 9, 12.5), y = MAX_VAL, label = lbl_abs,  parse = F, size=3) +  
            annotate("text", x = c(-12.5, -9, 0, 9, 12.5), y = MAX_VAL-(round((MAX_VAL * 5/100),1)), label = paste0("(", lbl_per, ")"),  parse = F, size=3) +
            theme_classic2() + scale_color_manual(values=c("#440154", "#21908c", "#fde725")) +
            theme(legend.position="top") + xlab("Distance to closest gene TSS in bp (log)") + ylab("Peak signalValue (Macs2; log2)") +
            labs(color = "Peak position")

# Marginal densities along x axis
xdens <- axis_canvas(p1, axis = "x") +
  geom_density(data = df, aes(x = dist_to_closest_gene, fill = annot.type),
               alpha = 0.7)+
  scale_fill_manual(values=c("#440154", "#21908c", "#fde725"))
# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
ydens <- axis_canvas(p1, axis = "y", coord_flip = TRUE)+
  geom_density(data = df, aes(x = signalValue, fill = annot.type),
               alpha = 0.7)+
  coord_flip()+
  scale_fill_manual(values=c("#440154", "#21908c", "#fde725"))
p1 <- insert_yaxis_grob(insert_xaxis_grob(p1, xdens, grid::unit(.2, "null"), position = "top"), ydens, grid::unit(.2, "null"), position = "right")
P <- annotate_figure(ggdraw(p1), top=title)

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## save output files
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
ggsave(filename = opt$outFile, plot = P, width = 9, height = 6)
q()
