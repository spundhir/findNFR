#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
	make_option(c("-i", "--inFile"), help="input file containing output from tss2binExpr script (can be stdin)"),
    make_option(c("-l", "--classVariable"), default="class", help="plot profiles categorized based on promoter a) class or b) directionality (default: %default)"),
    make_option(c("-c", "--directionalityCol"), default=8, help="column using which to define the promoter directionality classes (default: %default)"),
    make_option(c("-m", "--customGenes"), help="file containing custom class information (format: gene_name class) (default: all genes divided per directionality class"),
    make_option(c("-p", "--plotCol"), help="specific column for which profile needs to be plotted (default: all columns)"),
    make_option(c("-o", "--outFile"), help="output pdf file")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if(is.null(opt$inFile)) {
	cat("\nProgram: tss2binExpr2profileDirectionality.R (R script to plot profile directionality (eg. chd4))\n")
	cat("Author: BRIC, University of Copenhagen, Denmark\n")
	cat("Version: 1.0\n")
	cat("Contact: pundhir@binf.ku.dk\n");
	print_help(parser)
	q()
}

## load libraries
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(session))

## function to plot profile
plot_profile <- function(x, title, plot) {
    test.sub <- x

    if(missing(plot)) { plot="class"; }

    if(plot=="directionality") {
        test.melt <- melt(test.sub[,c(1:8,9,13)])

        col <- c(rep("#e31a1c", 3), rep("#1f78b4",2), rep("#33a02c",3))
        shp <- c(rep(17,3), rep(16,2), rep(18,3))

        ggplot(test.melt, aes(variable, group=1)) +
            stat_smooth(data = test.melt[which(test.melt$profile_directionality_class=="-1"),], aes(y = log(value), group=1, colour="-1"), method=lm, formula = y ~ poly(x,5), level=0.95) +
            stat_smooth(data = test.melt[which(test.melt$profile_directionality_class=="0"),], aes(y = log(value), group=1, colour="0"), method=lm, formula = y ~ poly(x,5), level=0.95) +
            stat_smooth(data = test.melt[which(test.melt$profile_directionality_class=="1"),], aes(y = log(value), group=1, colour="1"), method=lm, formula = y ~ poly(x,5), level=0.95) +
            ## colors are arranged alphabetically, meaning -1, 0 and 1
            scale_colour_manual("profile_directionality_class", breaks = c("-1", "0", "1"), values = c("#4faf4a","#ef7e1b", "#e31e1e")) +
            theme_bw() +
            facet_grid(.~class) +
            ggtitle(title) +
            scale_x_discrete(labels=c("-3", "-2", "-1", "HFR", "NFR", "1", "+2", "+3" )) +
            theme(axis.text.x = element_text(size=10, angle=45)) + xlab("") + ylab("Signal (TPKM, log)") +
            theme(legend.position="none")
    } else if(plot=="custom") {
        test.melt <- melt(test.sub[,c(1:8,15)])

        col <- c(rep("#e31a1c", 3), rep("#1f78b4",2), rep("#33a02c",3))   
        shp <- c(rep(17,3), rep(16,2), rep(18,3))        

        custom_class <- as.character(unique(test.melt$custom_class))

        ggplot(test.melt, aes(variable, group=1)) +              
            stat_smooth(data = test.melt[which(test.melt$custom_class==custom_class[1]),], aes(y = log(value), group=1, colour=custom_class[1]), method=lm, formula = y ~ poly(x,5), level=0.95) +              
            stat_smooth(data = test.melt[which(test.melt$custom_class==custom_class[2]),], aes(y = log(value), group=1, colour=custom_class[2]), method=lm, formula = y ~ poly(x,5), level=0.95) +
            stat_smooth(data = test.melt[which(test.melt$custom_class==custom_class[3]),], aes(y = log(value), group=1, colour=custom_class[3]), method=lm, formula = y ~ poly(x,5), level=0.95) +
            scale_colour_manual("custom_class", breaks = custom_class, values = c("#4faf4a","#ef7e1b", "#e31e1e")) +
            theme_bw() +
            ggtitle(title) +
            scale_x_discrete(labels=c("-3", "-2", "-1", "HFR", "NFR", "1", "2", "3" )) +
            theme(axis.text.x = element_text(size=10, angle=45)) + xlab("") + ylab("Signal (TPKM, log)") +
            theme(legend.position="top")
    } else {
        test.melt <- melt(test.sub[,c(1:8,9,11)])

        col <- c(rep("#e31a1c", 3), rep("#1f78b4",2), rep("#33a02c",3))
        shp <- c(rep(17,3), rep(16,2), rep(18,3))

        ggplot(test.melt, aes(variable, group=1)) +
            stat_smooth(data = test.melt[which(test.melt$promoter_class=="narrow"),], aes(y = log(value), group=1, colour="narrow"), method=lm, formula = y ~ poly(x,5), level=0.95) +
            stat_smooth(data = test.melt[which(test.melt$promoter_class=="medium"),], aes(y = log(value), group=1, colour="medium"), method=lm, formula = y ~ poly(x,5), level=0.95) +
            stat_smooth(data = test.melt[which(test.melt$promoter_class=="broad"),], aes(y = log(value), group=1, colour="broad"), method=lm, formula = y ~ poly(x,5), level=0.95) +
            ## colors are arranged alphabetically, meaning broad, medium and narrow
            scale_colour_manual("promoter_class", breaks = c("narrow", "medium", "broad"), values = c("#4faf4a","#ef7e1b", "#e31e1e")) +
            theme_bw() +
            facet_grid(.~class) +
            ggtitle(title) +
            scale_x_discrete(labels=c("-3", "-2", "-1", "HFR", "NFR", "1", "2", "3" )) +
            theme(axis.text.x = element_text(size=10, angle=45)) + xlab("") + ylab("Signal (TPKM, log)") +
            theme(legend.position="none")
    }
}


## read file
if(identical(opt$inFile, "stdin")==T) {
    data <- read.table(file("stdin"), comment.char = "")
} else {
    data <- read.table(opt$inFile, comment.char = "")
}

## assign column names
NCOL=ncol(data)
colnames(data) <- c("chr", "start", "end", "name", "score", "strand", "description")
colnames(data)[8:NCOL] <- seq(1,NCOL-7,1)

## format list containing profile values
class <- as.character(gsub("^.*#", "", data[grep("HFR", data$name),]$description))
distance <- as.numeric(apply(data[grep("HFR", data$name),], 1, function(x) unlist(strsplit(x[7], "#"))[7]))
gene <- apply(data[grep("HFR", data$name),], 1, function(x) unlist(strsplit(x[7], "#"))[4])
test <- lapply(data[,c(1:NCOL)], function(x) as.data.frame(matrix(x, nrow=length(x)/8, byrow = T)))
test <- lapply(test, setNames, nm <- c("2up", "1up", "up", "hfr", "nfr", "down", "1down", "2down"))
test <- lapply(test, cbind, class)
test <- lapply(test, cbind, distance)
test <- lapply(test, cbind, gene)

## promoter class
promoter_class <- test[[opt$directionalityCol]]$distance
promoter_class[which(test[[opt$directionalityCol]]$distance <= 500)] <- "narrow"
promoter_class[which(test[[opt$directionalityCol]]$distance > 500 & test[[opt$directionalityCol]]$distance <= 1000)] <- "medium"
promoter_class[which(test[[opt$directionalityCol]]$distance > 1000)] <- "broad"
test <- lapply(test, cbind, promoter_class)
data <- cbind(data[grep("HFR", data$name),], promoter_class)

## profile directionality
profile_directionality <- rep(0, nrow(test[[opt$directionalityCol]]))
profile_directionality <- (rowSums(test[[opt$directionalityCol]][c(5:8)])-rowSums(test[[opt$directionalityCol]][c(1:4)]))/(rowSums(test[[opt$directionalityCol]][c(5:8)])+rowSums(test[[opt$directionalityCol]][c(1:4)]))
test <- lapply(test, cbind, profile_directionality)
profile_directionality_class <- rep("0", nrow(test[[opt$directionalityCol]]))
profile_directionality_class[which(profile_directionality > 0.3)] <- "1"
profile_directionality_class[which(profile_directionality < -0.3)] <- "-1"
test <- lapply(test, cbind, profile_directionality_class)
data <- cbind(data[grep("HFR", data$name),], profile_directionality)
data <- cbind(data[grep("HFR", data$name),], profile_directionality_class)
data$name <- gsub("HFR#", "", data$name)

## filter based on input genes and define custom classes if provided
if(!is.null(opt$customGenes)) {
    customGenes <- read.table(opt$customGenes)
    colnames(customGenes) <- c("gene", "custom_class")

    # slow but genes in the same order as in input (opt$customGenes)
    #lapply(test, function(x) x[unlist(lapply(as.character(customGenes$gene), function(y) grep(sprintf("^%s$", y),x$gene))),])
    # fast but genes not in the same order as in input (opt$customGenes)
    test <- lapply(test, function(x) x[grep(paste(sprintf("^%s$", as.character(customGenes$gene)), collapse="|"), x$gene),])
    test <- lapply(test, function(x) { x <- merge(x, customGenes, by.x="gene", by.y="gene")[,c(2:ncol(x),1,ncol(x)+1)]; colnames(x)[ncol(x)] <- "custom_class"; return(x); })
    data <- data[grep(paste(sprintf("^%s$", as.character(customGenes$gene)), collapse="|"), data$name),]
}

## plot profiles or output file to stdout
if(!is.null(opt$outFile)) {
    if(!is.null(opt$plotCol)) {
        factors_order <- names(test)[as.numeric(opt$plotCol)]
    } else {
        factors_order <- names(test)[c(8:NCOL)]
    }

    if(!is.null(opt$customGenes)) {
        p <- lapply(factors_order, function(x) { plot_profile(test[[x]], x, "custom") })
    } else {
        p <- lapply(factors_order, function(x) { plot_profile(test[[x]], x, opt$classVariable) })
    }

    g <- do.call(grid.arrange, p)
    ggsave(g, dpi=300, height=30, width=25, filename=opt$outFile, useDingbats=FALSE)
} else {
    write.table(data, "", sep="\t", col.names=T, row.names=F, quote=F)
}

## write output session file
#outSession <- sprintf("%s.session", opt$outFile)
#save.session(outSession)
q()

