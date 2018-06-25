#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
	make_option(c("-i", "--inFile"), help="input file containing output from tss2binExpr script (can be stdin)"),
    make_option(c("-c", "--profileCol"), default=8, help="column using which to define the profile directionality classes (default: %default)"),
    make_option(c("-o", "--outFile"), help="output pdf file"),
    make_option(c("-l", "--outFileClass"), default="promoter_class", help="plot profiles categorized based on \"promoter_class\" or \"profile_directionality_class\" (default: %default)")
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

    if(missing(plot)) { plot="promoter_class"; }

    if(plot=="profile_directionality_class") {
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
test <- lapply(data[,c(1:NCOL)], function(x) as.data.frame(matrix(x, nrow=length(x)/8, byrow = T)))
test <- lapply(test, setNames, nm <- c("2up", "1up", "up", "hfr", "nfr", "down", "1down", "2down"))
test <- lapply(test, cbind, class)
test <- lapply(test, cbind, distance)

## promoter class
promoter_class <- test[[opt$profileCol]]$distance
promoter_class[which(test[[opt$profileCol]]$distance <= 500)] <- "narrow"
promoter_class[which(test[[opt$profileCol]]$distance > 500 & test[[opt$profileCol]]$distance <= 1000)] <- "medium"
promoter_class[which(test[[opt$profileCol]]$distance > 1000)] <- "broad"
test <- lapply(test, cbind, promoter_class)
data <- cbind(data[grep("HFR", data$name),], promoter_class)

## profile directionality
profile_directionality <- rep(0, nrow(test[[opt$profileCol]]))
profile_directionality <- (rowSums(test[[opt$profileCol]][c(5:8)])-rowSums(test[[opt$profileCol]][c(1:4)]))/(rowSums(test[[opt$profileCol]][c(5:8)])+rowSums(test[[opt$profileCol]][c(1:4)]))
test <- lapply(test, cbind, profile_directionality)
profile_directionality_class <- rep("0", nrow(test[[opt$profileCol]]))
profile_directionality_class[which(profile_directionality > 0.3)] <- "1"
profile_directionality_class[which(profile_directionality < -0.3)] <- "-1"
test <- lapply(test, cbind, profile_directionality_class)
data <- cbind(data[grep("HFR", data$name),], profile_directionality)
data <- cbind(data[grep("HFR", data$name),], profile_directionality_class)

## plot profiles
if(!is.null(opt$outFile)) {
    factors_order <- names(test)[c(8:NCOL)]
    p <- lapply(factors_order, function(x) { plot_profile(test[[x]], x, opt$outFileClass) })
    g <- do.call(grid.arrange, p)
    ggsave(g, dpi=300, height=30, width=25, filename=sprintf("%s.pdf", opt$outFile), useDingbats=FALSE)
}

## write output file
data$name <- gsub("HFR#", "", data$name)
write.table(data, "", sep="\t", col.names=T, row.names=F, quote=F)

## write output session file
#outSession <- sprintf("%s.session", opt$outFile)
#save.session(outSession)
#q()

