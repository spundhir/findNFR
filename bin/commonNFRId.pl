#!/usr/bin/perl -w

=copyright_info
commonNFRId.pl: determine common Nucleosome Free Regions (NFR) between two replicates
Copyright (C) 2015  Sachin Pundhir (pundhir@binf.ku.dk)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
=cut

use strict;
use warnings;
use Getopt::Long;
use Tie::IxHash;
#use List::MoreUtils qw(uniq);

###############################################################################
## parse input options
use vars qw($nfrFiles $bamFiles $sizeFactors $extends $genome $help);

$genome="mm9";

GetOptions ("i=s"  => \$nfrFiles,
            "k=s"  => \$bamFiles,
            "m=s"  => \$sizeFactors,
            "c=s"  => \$extends,
            "y=s"  => \$genome,
            "help" => \$help,
            "h"    => \$help);

usage() if($help || !$nfrFiles || !$bamFiles || !$sizeFactors);

###############################################################################
sub usage {
	print STDERR "\nProgram: commonNFRId.pl (determine common Nucleosome Free Regions (NFR) between multiple replicates based on identifier)\n";
	print STDERR "Author: BRIC, University of Copenhagen, Denmark\n";
	print STDERR "Version: 1.0\n";
	print STDERR "Contact: pundhir\@binf.ku.dk\n";
	print STDERR "Usage: commonNFR.pl -i <file>\n";
	print STDERR " -i <file>         [input file(s) having NFRs]\n";
    print STDERR "                   [if multiple, please separate them by a comma]\n";
    print STDERR " -k <file>         [input BAM file(s)]\n";
    print STDERR "                   [if multiple, please separate them by a comma]\n";
    print STDERR " -m <float>        [size factor(s) to normalize expression of reads]\n";
    print STDERR "                   [if multiple, please separate them by a comma]\n";
    print STDERR "[OPTIONS]\n";
    print STDERR " -c <int>          [value(s) to extend 3' end of reads by input number of bases (default: 0)]\n";
    print STDERR "                   [if multiple, please separate them by a comma]\n";
    print STDERR " -y <string>       [genome (default: mm9)]\n";
	print STDERR " -h                [help]\n\n";
	exit(-1);
}

###############################################################################

## populate genome file based on input genome
my $GENOME_FILE=`initialize_genome -i $ENV{FINDNFRPATH}/data/annotations/GENOME_FILE -g $genome`;
#$GENOME_FILE="$ENV{FINDNFRPATH}/data/annotations/$GENOME_FILE";
if(! -f $GENOME_FILE) {
    print "\ncomputation for $genome is not feasible yet\n";
    print "please add the chromosome size file for $genome at $ENV{FINDNFRPATH}/data/annotations\n";
    print "also update the $ENV{FINDNFRPATH}/data/annotations/GENOME_FILE\n";
    usage();
}

## set default value to extends parameter, if not defined
my @nfrFileArray=split(/\,/,$nfrFiles);
if(!defined($extends)) {
    $extends="";
    foreach(@nfrFileArray) {
        $extends.="0,";
    }
    $extends=~s/\,$//g;
}

$nfrFiles=~s/\,/ /g;
chomp($nfrFiles);

my @data=();
@data = `zless $nfrFiles | sort -k 4,4 -k 15r,15`;

## determine unqiue highest scoring (significant) NFRs from multiple replicates
my %NFR=(); my @F=();
tie %NFR, 'Tie::IxHash';
foreach(@data) {
    chomp($_);
    @F=split(/\s+/,$_);
    if($F[14]=~/SIG_PEAKS/) { $F[14]=1; }
    else { $F[14]=0; }

    if(!defined($NFR{$F[3]})) {
        $NFR{$F[3]}{'info'}=$_;
        $NFR{$F[3]}{'count'}=1;
        $NFR{$F[3]}{'score'}=$F[4];
        $NFR{$F[3]}{'significant'}=$F[14];
    } elsif($F[4]!~/NA/ && $F[14] > $NFR{$F[3]}{'significant'}) {
        $NFR{$F[3]}{'info'}=$_;
        $NFR{$F[3]}{'count'}++;
        $NFR{$F[3]}{'score'}=$F[4];
        $NFR{$F[3]}{'significant'}=$F[14];
    } elsif($F[4]!~/NA/ && $F[14] >= $NFR{$F[3]}{'significant'} && $F[4] > $NFR{$F[3]}{'score'}) {
        $NFR{$F[3]}{'info'}=$_;
        $NFR{$F[3]}{'count'}++;
        $NFR{$F[3]}{'score'}=$F[4];
        $NFR{$F[3]}{'significant'}=$F[14];
    }
    #push(@{$NFR{$F[3]}{'files'}}, $F[15]);
}

## recompute nfrDip score using all input bam files
foreach(keys(%NFR)) {
    chomp($NFR{$_}{'info'});
    @F=split(/\s+/,$NFR{$_}{'info'});

    if($F[1]>0 && $F[2]>0 && $F[4]!~/NA/) { 
        $F[8]=`coor2expr -i $F[6] -j $bamFiles -k $sizeFactors -d -e $extends -g $genome`;
        $F[12]=`coor2expr -i $F[0]:$F[1]-$F[2] -j $bamFiles -k $sizeFactors -d -e $extends -g $genome`;
        $F[10]=`coor2expr -i $F[7] -j $bamFiles -k $sizeFactors -d -e $extends -g $genome`;
        $F[4]=((($F[8]+$F[10])/($F[9]+$F[11]))-($F[12]/$F[13]));
        $F[4]=sprintf("%0.4f", $F[4]);
    }
    my $line="";
    foreach(@F) {
        chomp($_);
        $line.="$_\t";
    }
    $line=~s/\t$//g;
    print "$line\n";
    #print "$NFR{$_}{'info'}\t". uniq(@{$NFR{$_}{'files'}})."\n";
}

exit;
