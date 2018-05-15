#!/usr/bin/env perl
use strict;
use Math::Round;

## By Andrew Spriggs, CSIRO Ag&Food, 2018 ##
## andrew.spriggs@csiro.au ##
## https://github.com/spriggsy83 ##

if($#ARGV < 0){
	warn("Produces alignment statistics from log-file output of 'biokanga align'\n");
	die("Usage: $0 <in align log file> [ more align log files ]\n");
}

my @totReads;
my @uniqueHits;
my @multiHits;
my @noHits;
my @inFiles;

for(my $i=0; $i<=$#ARGV; $i++){
	push @inFiles, $ARGV[$i];
	push @totReads, 0;
	push @uniqueHits, 0;
	push @multiHits, 0;
	push @noHits, 0;
}

for(my $fileNum = 0; $fileNum <= $#inFiles; $fileNum++){
	open(INFILE, "<".$inFiles[$fileNum]) or die("Couldn't open ".$inFiles[$fileNum]." for input\n");
	
	while(<INFILE>){
		chomp;
		if(/From (\d+) source reads there are/){
			$totReads[$fileNum] = $1;
		}elsif(/(\d+) \(AA\) Alignment accepted/){
			$uniqueHits[$fileNum] = $1;
		}elsif(/(\d+) \(ML\) Aligned to multiloci/){
			$multiHits[$fileNum] = $1;
		}elsif(/(\d+) \(NL\) No potential alignment loci/){
			$noHits[$fileNum] = $1;
		}
	}	
	
	close INFILE;
}

print "\tTotalReads\tUniquelyAligned\tMultiLoci\tNoAlignLoci\tOtherNoRes\tUniquelyAligned\tMultiLoci\tNoAlignment\n";
for(my $fileNum = 0; $fileNum <= $#inFiles; $fileNum++){
	print $inFiles[$fileNum]."\t";
	print $totReads[$fileNum]."\t";
	print $uniqueHits[$fileNum]."\t";
	print $multiHits[$fileNum]."\t";
	print $noHits[$fileNum]."\t";
	my $remain = $totReads[$fileNum] - $uniqueHits[$fileNum] - $multiHits[$fileNum] - $noHits[$fileNum];
	print $remain."\t";
	if($totReads[$fileNum] != 0){
		print nearest(0.1, $uniqueHits[$fileNum] / $totReads[$fileNum] * 100)."%\t";
		print nearest(0.1, $multiHits[$fileNum] / $totReads[$fileNum] * 100)."%\t";
		print nearest(0.1, ($noHits[$fileNum] + $remain) / $totReads[$fileNum] * 100)."%\n";
	}else{
		print "0%\t0%\t0%\n";
	}
}
print "\n\n";

