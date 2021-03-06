#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;

###############################################################################
## parse input options
use vars qw($inFile $db $help);
$db="NA";

GetOptions ("i=s"  => \$inFile,
            "j=s"  => \$db,
            "help" => \$help,
            "h"    => \$help);

usage() if($help);

###############################################################################
sub usage {
	print STDERR "\nProgram: meme2homer.pl (convert from MEME to HOMER format)\n";
	print STDERR "Author: BRIC, University of Copenhagen, Denmark\n";
	print STDERR "Version: 1.0\n";
	print STDERR "Contact: pundhir\@binf.ku.dk\n";
	print STDERR "Usage: meme2homer.pl -i <file> [OPTIONS]\n";
	print STDERR " -i <file>         [input file in MEME format (can be stdin)]\n";
    print STDERR "[OPTIONS]\n";
    print STDERR " -j <string>       [database name (default: NA)]\n";
	print STDERR " -h                [help]\n";
    print STDERR "[NOTE]\n";
    print STDERR " Inspired from: http://homer.ucsd.edu/homer/motif/creatingCustomMotifs.html\n";
	exit(-1);
}

###############################################################################

my @data=();

if(!defined($inFile)) {
    usage();
}
elsif($inFile=~/^stdin$/) {
    my $INFILE=();
    $INFILE=*STDIN;
    @data=<$INFILE>;
}
else {
    chomp($inFile);
    @data=`zless $inFile`;
}

my $id=(); my $des=();
my $start=(); my $score=();
my @freq=(); my @matrix=();
my $max=(); my @F=();
my @values=();

foreach my $l(@data) {
    @F=split(/\s+/,$l);
    if($l=~/MOTIF/) {
        if(defined($F[2])) {
            $id=$F[1]; $des="$F[2]/$db";
        }
        else {
            $id=$F[1]; $des="$F[1]/$db";
        }
    }
    elsif($l=~/letter-pro/) {
        $start=1; $score=0;
        @freq=(); @matrix=();
    }
    elsif($start && ($l=~/^$/ || $l=~/[a-zA-Z]+/)) {
        foreach(@freq) {
            $score+=log($_/0.25);
        }
        $score=$score-4;
        print ">$id\t$des\t$score\n";
        foreach(@matrix) {
            @values=split(/\s+/,$_);
            print "$values[0]\t$values[1]\t$values[2]\t$values[3]\n";
        }
        $start=0;
    }
    elsif($start) {
        $l=~s/^\s+//g;
        push(@matrix, $l);
        $max=0;
        my @T=split(/\s+/,$l);
        foreach(@T) {
            if($_>$max) {
                $max=$_;
            }
        }
        push(@freq, $max);
    }
}

## print output for the last motif read
if($start) {
    foreach(@freq) {
        $score+=log($_/0.25);
    }
    $score=$score-4;
    print ">$id\t$des\t$score\n";
    foreach(@matrix) {
        @values=split(/\s+/,$_);
        print "$values[0]\t$values[1]\t$values[2]\t$values[3]\n";
    }
    $start=0;
}

exit;
