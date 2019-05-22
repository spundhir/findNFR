#!/usr/bin/perl -w

use strict;
use warnings;
use Tie::IxHash;
use Getopt::Long;
use Statistics::Basic qw(:all);
use perlModule;

###############################################################################
## parse input options
use vars qw($regionFile $bamFile $outDir $sizeFactor $pvalue $distribution $minClusterHeight $minBlockHeight $distance $scale $blockHeight $minNFRLength $maxNFRLength $extend $genome $fileSuffix $onlyUnique $help);
$pvalue=0.05;
$distribution="pois";
$distance=70;
$scale=0.6;
$minClusterHeight=2;
$minBlockHeight=2;
$blockHeight="abs";
$minNFRLength=20;
$maxNFRLength=1000;
$extend=0;
$genome="mm9";
$fileSuffix="";

GetOptions ("s=s"  => \$regionFile,
            "b=s"  => \$bamFile,
            "o=s"  => \$outDir,
            "z=s"  => \$sizeFactor,
            "p=s"  => \$pvalue,
            "d=s"  => \$distribution,
            "x=s"  => \$distance,
            "l=s"  => \$scale,
            "g=s"  => \$blockHeight,
            "n=s"  => \$minNFRLength,
            "v=s"  => \$maxNFRLength,
            "e=s"  => \$extend,
            "y=s"  => \$genome,
            "f=s"  => \$fileSuffix,
            "u"    => \$onlyUnique,
            "help" => \$help,
            "h"    => \$help);

usage() if($help || !$regionFile || !$bamFile || !$outDir || !$sizeFactor);

###############################################################################
sub usage {
	print STDERR "\nProgram: findNFR.pl (determine Nucleosome Free Regions (NFR) using ChIP-seq data for histone marks)\n";
	print STDERR "Author: BRIC, University of Copenhagen, Denmark\n";
	print STDERR "Version: 1.0\n";
	print STDERR "Contact: pundhir\@binf.ku.dk\n";
	print STDERR "Usage: findNFR.pl -s <file> -b <file> -o <dir> -z <float> [OPTIONS]\n";
	print STDERR " -s <file>         [file containing regions of interest in BED format]\n";
    print STDERR "                   [if multiple, seperate them by a comma]\n";
	print STDERR " -b <file>         [histone ChIP-seq file in BAM format]\n";
	print STDERR " -o <dir>          [directory where output files will be kept]\n";
	print STDERR " -z <float>        [size factor to normalize read expression]\n";
	print STDERR "[OPTIONS]\n";
    print STDERR " -p <float>        [pvalue for read enrichment (default: 0.05)]\n";
    print STDERR " -d <string>       [distribution to use for data fit (pois or nbinom) (default: pois)]\n";
	print STDERR " -x <int>          [maximum distance between the blocks (default: 70)]\n";
	print STDERR " -l <float>        [scale to define blocks (default: 0.6)]\n";
	print STDERR " -g <string>       [relative block height (abs or rel) (default: abs)]\n";
    print STDERR " -n <int>          [minimum length of nucleosome free region (default: 20)]\n";
    print STDERR " -v <int>          [maximum length of nucleosome free region (default: 1000)]\n";
    print STDERR " -e <int>          [extend 3' end of reads by input number of bases (default: 0)]\n";
    print STDERR " -y <string>       [genome (default: mm9)]\n";
    print STDERR " -f <string>       [a string added at the end of output files. useful when running in parallel]\n";
    print STDERR " -u                [output one NFR having highest nfrDip score corresponding to each region]\n";
	print STDERR " -h                [help]\n\n";
	exit(-1);
}

###############################################################################

my $ID=$bamFile;
$ID=~s/^.*\///g;
$ID=~s/\.gz$//g;
my $start=(); my $end=(); my $coor=(); my @data=(); my $COUNT=0;

## create output directory, if does not exist
if ( ! -d $outDir) {
    system("mkdir $outDir");
}

## estimate minimum block and block group height, define blocks and block groups in histone enriched (summit) region and define NFRs
$regionFile=~s/\,/ /g;
chomp($regionFile);
@data=`zless $regionFile | sortBed -i stdin`;

# remove nfr file, if already exists
if(-e "$outDir/$ID.nfr$fileSuffix") {
    system("rm $outDir/$ID.nfr$fileSuffix");
}

foreach my $l(@data) {
    my @F=split(/\s+/,$l);
    my $coor="$F[0]:$F[1]-$F[2]";
    ## Step-1: retrieve reads corresponding to input coordinate
    system("samtools view -b $bamFile $coor | bedtools bamtobed -i - | perl -ane '\$F[5]=~s/\-\$/+/g; \$F[4]=sprintf(\"%0.2f\", 1/$sizeFactor); print \"\$F[0]\\t\$F[1]\\t\$F[2]\\t\$F[3]\\t\$F[4]\\t\$F[5]\\n\";' > $outDir/$ID.tmp$fileSuffix");

    ## break, if minimum 10 reads corresponding to the coordinate are not available
    $COUNT=`cat $outDir/$ID.tmp$fileSuffix | wc -l`;

    if($COUNT>=10) {
        ## Step-2: define block group and blocks corresponding to input coordinate
        $distance=$F[2]-$F[1]; # define only one block group for each region
        system("reads2peaks -i $outDir/$ID.tmp$fileSuffix -o optimizeThreshold -p $pvalue -d $distribution -x $distance -w $minNFRLength > $outDir/$ID.bg$fileSuffix");

        ## break, if minimum 2 significant peaks (block groups) corresponding to the coordinate are not available
        my $COUNT=`cat $outDir/$ID.bg$fileSuffix | wc -l`;

        ## Step-3: define nuclesome free regions
        if($COUNT>=2) {
            if(defined($onlyUnique)) {
                system("peaks2NFR.pl -i $outDir/$ID.bg$fileSuffix -b $bamFile -z $sizeFactor -n $minNFRLength -v $maxNFRLength -e $extend -y $genome | perl -ane 'chomp(\$_); print \"\$_\\tSIG_PEAKS\\n\";' | perl -ane 'BEGIN { \$max=0; } if(\$F[4]>\$max) { chomp(\$_); \$max_nfr=\$_; \$max=\$F[4]; } END { print \"\$max_nfr\\n\"; }' >> $outDir/$ID.nfr$fileSuffix");
            }
            else {
                system("peaks2NFR.pl -i $outDir/$ID.bg$fileSuffix -b $bamFile -z $sizeFactor -n $minNFRLength -v $maxNFRLength -e $extend -y $genome | perl -ane 'chomp(\$_); print \"\$_\\tSIG_PEAKS\\n\";' >> $outDir/$ID.nfr$fileSuffix");
            }
        }
        else {
            print STDERR "No significant peak (block group) is found for $coor, instead computing NFR(s) for all blocks\n";
            system("reads2peaks -i $outDir/$ID.tmp$fileSuffix -o optimizeThreshold -p $pvalue -d $distribution -x $distance -w $minNFRLength -a > $outDir/$ID.bg$fileSuffix");

            ## break, if minimum 2 peaks (block groups) corresponding to the coordinate are not available
            my $COUNT=`cat $outDir/$ID.bg$fileSuffix | wc -l`;

            if($COUNT>=2) { 
                if(defined($onlyUnique)) {
                    system("peaks2NFR.pl -i $outDir/$ID.bg$fileSuffix -b $bamFile -z $sizeFactor -n $minNFRLength -v $maxNFRLength -e $extend -y $genome | perl -ane 'chomp(\$_); print \"\$_\\tNONSIG_PEAKS\\n\";' | perl -ane 'BEGIN { \$max=0; } if(\$F[4]>\$max) { chomp(\$_); \$max_nfr=\$_; \$max=\$F[4]; } END { print \"\$max_nfr\\n\"; }' >> $outDir/$ID.nfr$fileSuffix");
                }
                else {
                    system("peaks2NFR.pl -i $outDir/$ID.bg$fileSuffix -b $bamFile -z $sizeFactor -n $minNFRLength -v $maxNFRLength -e $extend -y $genome | perl -ane 'chomp(\$_); print \"\$_\\tNONSIG_PEAKS\\n\";' >> $outDir/$ID.nfr$fileSuffix");
                }
            }
            else {
                print STDERR "No peak (block group) is found for $coor, putting dummy numbers\n";
                system("echo -e $F[0]\t$F[1]\t$F[2]\t$coor\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA >> $outDir/$ID.nfr$fileSuffix");
            }
        }
    }
    else {
        print STDERR "No reads found for $coor, putting dummy numbers\n";
        system("echo -e $F[0]\t$F[1]\t$F[2]\t$coor\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA >> $outDir/$ID.nfr$fileSuffix");

    }
}
system("rm $outDir/$ID.tmp$fileSuffix");
system("rm $outDir/$ID.bg$fileSuffix");

## print output file created above
print "$outDir/$ID.nfr$fileSuffix\n";

exit(0);
