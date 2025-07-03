#!/usr/bin/env Rscript
#suppressPackageStartupMessages(library("optparse", lib.loc = "/home/users/sachin/R/x86_64-redhat-linux-gnu-library/2.15"))
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
  make_option(c("-i", "--inDistFile"), help="input file containing peak score and distance to TSS information"),
  make_option(c("-j", "--inJasparFile"), help="input file containing motifs enriched at peaks"),
  make_option(c("-o", "--outFile"), help="output pdf file")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if(((is.null(opt$inDistFile) & is.null(opt$inJasparFile)) | is.null(opt$outFile))) {
  cat("\nProgram: peak2qc.R (plot peak QC information)\n")
  cat("Author: BRIC, University of Copenhagen, Denmark\n")
  cat("Version: 1.0\n")
  cat("Contact: pundhir@binf.ku.dk\n");
  print_help(parser)
  q()
}

suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("ggrastr"))
suppressPackageStartupMessages(library("ggpubr"))

if(!is.null(opt$inDistFile)) {
    ## read input file
    if(identical(opt$inDistFile, "stdin")==T) {
        df <- read.table(file("stdin"))
    } else {
        df <- read.table(opt$inDistFile)
    }
    colnames(df) <- c("chr", "start", "end", "name", "score", "strand", "signalValue", "dist_to_geneTSS")

    ## reorganize values to plot
    df$dist_to_geneTSS <- log(abs(df$dist_to_geneTSS)+1) * ifelse(df$dist_to_geneTSS<0, -1, 1)
    df$annot.type <- "TSS"
    df[which(abs(df$dist_to_geneTSS) > log(1) & abs(df$dist_to_geneTSS) <= log(50000)),]$annot.type <- "proximal"
    df[which(abs(df$dist_to_geneTSS) > log(50000)),]$annot.type <- "distal"
    df$annot.type <- factor(df$annot.type, levels=c("TSS", "proximal", "distal"), ordered=T)
    table(df$annot.type)

    brk <- c(min(df$dist_to_geneTSS), -log(50000), 0, log(50000), max(df$dist_to_geneTSS))
    lbl <- paste(sprintf("%0.1f", (table(cut(df$dist_to_geneTSS, 
                breaks=c(min(df$dist_to_geneTSS), -log(50000), -log(1.01), 0, log(1.01), log(50000), max(df$dist_to_geneTSS)), 
                labels=c(min(df$dist_to_geneTSS), -log(50000), 0, 0, log(50000), max(df$dist_to_geneTSS))))*100)/nrow(df)), "%")

    ## make the plot
    if(length(which(!is.na(df$signalValue))) > 0) {
        p1 <- ggplot(df, aes(dist_to_geneTSS, signalValue)) 
    } else {
        p1 <- ggplot(df, aes(dist_to_geneTSS, score)) 
    }

    p1 <- p1 + geom_point_rast(aes(color=annot.type), alpha=0.2) +
                stat_density_2d(geom = "polygon", contour = TRUE, aes(fill = after_stat(level)), colour = "black", alpha=0.5, show.legend = F) +
                scale_fill_distiller(palette = "Reds", direction = 1) + 
                #xlim(c(-15,15)) + 
                geom_vline(xintercept = brk, lty=2) +
                scale_x_continuous(breaks=brk, labels = as.numeric(sprintf("%0.0f", brk))) +
                #annotate("text", x = c(-10, -5, 0, 5, 10), y = 2, label = lbl,  parse = F, size=3) +  theme_classic2() +
                annotate("text", x = c(-12.5, -9, 0, 9, 12.5), y = 2, label = lbl,  parse = F, size=3) +  theme_classic2()
                theme(legend.position="top") + xlab("Distance to closest gene TSS in bp (log)") + ylab("Peak signalValue (Macs2)") +
                labs(color = "Peak position")
} else if(!is.null(opt$inJasparFile)) {
    df <- read.table(opt$inJasparFile)
}
ggsave(filename = opt$outFile, plot = p, width = 4, height = 3)
q()
