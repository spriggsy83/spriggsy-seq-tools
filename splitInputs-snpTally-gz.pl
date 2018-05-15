#!/usr/bin/env perl
use IO::Compress::Gzip qw(gzip $GzipError) ;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
use strict;

## By Andrew Spriggs, CSIRO Ag&Food, 2018 ##
## andrew.spriggs@csiro.au ##
## https://github.com/spriggsy83 ##

if($#ARGV != 3){
	warn "Usage: $0 <in seqLengths> <in refseqs> <in alignList file> <num splits>\n";
	warn "AlignList file is tab-seprated:\n";
	warn "\tSampleName\tSAMfile\tSNPfile\n";
	warn "Fasta, SAM and SNP table files expected as .gz compressed\n";
	die("See README-tallySNPs.txt\n");
}

my $inSeqlensFile = $ARGV[0];
my $inSeqsFile = $ARGV[1];
my $inListFile = $ARGV[2];
my $numSplits = $ARGV[3];

if(!($inSeqsFile =~ /\.gz$/)){
	die("Fasta, SAM and SNP table files expected as .gz compressed\n");
}

my %refSeqs;
my @inIDs;
my @inSAMFiles;
my @inSNPFiles;

open(INFILE, "<$inListFile") or die("Couldn't open $inListFile for input\n");
while(<INFILE>){
	chomp;
	my (@line) = (split /\t/, $_);
	if($#line == 2){
		if(!($line[1] =~ /\.gz$/ && $line[2] =~ /\.gz$/)){
			die("Fasta, SAM and SNP table files expected as .gz compressed\n");
		}
		push @inIDs, $line[0];
		push @inSAMFiles, $line[1];
		push @inSNPFiles, $line[2];
	}
}
close INFILE;

open(INFILE, "<$inSeqlensFile") or die("Couldn't open $inSeqlensFile for input\n");
while(<INFILE>){
	chomp;
	my (@line) = (split /\t/, $_);
	if($#line >= 1){
		$refSeqs{$line[0]} = 1;
	}
}
close INFILE;

my $totalRefSeqs = scalar keys %refSeqs;
while($totalRefSeqs % $numSplits != 0){
	$totalRefSeqs++;
}
my $numPerFile = $totalRefSeqs / $numSplits;
my %refFileNums;
my $count = 0;
my $fileI = 0;
foreach my $refID (keys %refSeqs){
	$refFileNums{$refID} = $fileI;
	$count++;
	if($count == $numPerFile){
		$count = 0;
		$fileI++;
	}
}

`mkdir -p splitInputs`;
`mkdir -p splitInputs/sam`;
`mkdir -p splitInputs/snp`;
`mkdir -p splitInputs/fasta`;
`mkdir -p splitInputs/fLists`;
`mkdir -p splitInputs/sLists`;


	## File list split ## 
for(my $i=0; $i < $numSplits; $i++){
	open(OUTFILE, ">splitInputs/fLists/inFilesList.pt$i.txt") or die("Couldn't open splitInputs/fLists/inFilesList.pt$i.txt for output\n");
	for(my $y=0; $y <= $#inIDs; $y++){
		print OUTFILE $inIDs[$y]."\tsplitInputs/sam/".$inIDs[$y].".pt$i.sam.gz\tsplitInputs/snp/".$inIDs[$y].".pt$i.snps.gz\n";
	}
	close OUTFILE;
}
	

	## Ref list files split ## 
my @outSLenFiles;
for(my $i=0; $i < $numSplits; $i++){
	open my $file, ">splitInputs/sLists/refSeqList.pt$i.txt" or die("Couldn't open splitInputs/sLists/refSeqList.pt$i.txt for output\n");
	push @outSLenFiles, $file;
}

open(INFILE, "<$inSeqlensFile") or die("Couldn't open $inSeqlensFile for input\n");
print "Splitting $inSeqlensFile...\n";
while(<INFILE>){
	chomp;
	my (@line) = (split /\t/, $_);
	if($#line >= 1){
		print { $outSLenFiles[ $refFileNums{$line[0]} ] } $_."\n";
	}
}
close INFILE;

foreach my $lensFile (@outSLenFiles){
	close $lensFile;
}


	## Ref seq files split ## 
my @outSeqsFiles;
for(my $i=0; $i < $numSplits; $i++){
	my $outfile = new IO::Compress::Gzip "splitInputs/fasta/refSeqs.pt$i.fasta.gz"
		or die "Couldn't open splitInputs/fasta/refSeqs.pt$i.fasta.gz for output.gz\nIO::Compress::Gzip failed: $GzipError\n";
	push @outSeqsFiles, $outfile;
}

my $INFILE = new IO::Uncompress::Gunzip $inSeqsFile
	or die "Couldn't open $inSeqsFile for input.\nIO::Uncompress::Gunzip failed: $GunzipError\n";
print "Splitting $inSeqsFile...\n";

my $refID = "";
while(<$INFILE>){
	chomp;
	if(/^>(\S+)/){
		$refID = $1;
		if(defined($refFileNums{$refID})){
			print { $outSeqsFiles[ $refFileNums{$refID} ] } ">$refID\n";
		}else{
			$refID = "";
		}
	}else{
		if($refID ne ""){
			print { $outSeqsFiles[ $refFileNums{$refID} ] } $_."\n";
		}
	}
}
close $INFILE;

foreach my $seqsFile (@outSeqsFiles){
	close $seqsFile;
}


	# SNP files split ## 
for(my $y=0; $y <= $#inIDs; $y++){
	my @outSnpFiles;
	for(my $i=0; $i < $numSplits; $i++){
		my $snpFile = new IO::Compress::Gzip "splitInputs/snp/".$inIDs[$y].".pt$i.snps.gz"
			or die "Couldn't open splitInputs/snp/".$inIDs[$y].".pt$i.snps.gz for output.gz\nIO::Compress::Gzip failed: $GzipError\n";
		push @outSnpFiles, $snpFile;
	}
	
	my $INFILE = new IO::Uncompress::Gunzip $inSNPFiles[$y]
		or die "Couldn't open ".$inSNPFiles[$y]." for input.\nIO::Uncompress::Gunzip failed: $GunzipError\n";
	print "Splitting ".$inSNPFiles[$y]."...\n";
	while(<$INFILE>){
		chomp;
		my (@line) = (split /,/, $_);
		if($#line > 4){
			$line[3] =~ s/"//g;
			if( $line[3] ne "Chrom" ){
				if(defined($refFileNums{$line[3]})){
					print { $outSnpFiles[ $refFileNums{$line[3]} ] } $_."\n";
				}
			}
		}
	}
	close $INFILE;
	
	foreach my $snpFile (@outSnpFiles){
		close $snpFile;
	}
}


	## SAM files split ## 
for(my $y=0; $y <= $#inIDs; $y++){
	my @outSamFiles;
	for(my $i=0; $i < $numSplits; $i++){
		my $samFile = new IO::Compress::Gzip "splitInputs/sam/".$inIDs[$y].".pt$i.sam.gz"
			or die "Couldn't open splitInputs/sam/".$inIDs[$y].".pt$i.sam.gz for output.gz\nIO::Compress::Gzip failed: $GzipError\n";
		push @outSamFiles, $samFile;
	}
	
	my $INFILE = new IO::Uncompress::Gunzip $inSAMFiles[$y]
		or die "Couldn't open ".$inSAMFiles[$y]." for input.\nIO::Uncompress::Gunzip failed: $GunzipError\n";
	print "Splitting ".$inSAMFiles[$y]."...\n";
	while(<$INFILE>){
		chomp;
		my (@line) = (split /\t/, $_);
		if($#line > 4){
			if(defined($refFileNums{$line[2]})){
				print { $outSamFiles[ $refFileNums{$line[2]} ] } $_."\n";
			}
		}
	}
	close $INFILE;
	
	foreach my $samFile (@outSamFiles){
		close $samFile;
	}
}


