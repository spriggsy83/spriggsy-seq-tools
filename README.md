## README
This is a collection of small C++ programs and Perl scripts that have been handy in bioinformatic sequence analysis.
Most parse/manipulate fasta/fastq formatted DNA/RNA sequence files.  Some read SAM format sequence alignment results.
A couple deal specifically with results from the program 'biokanga align' (https://github.com/csiro-crop-informatics/biokanga).

C++ program prerequisites: Boost C++ Libraries (www.boost.org) and OpenMPI (www.open-mpi.org/)

Compile C++ programs on Unix systems with:
`bash buildAllC.sh`

|Program|Use|
|--|--|
|getSubSeqs|Extract sub-sequences, given a list of coordinate ranges|
|extractSeqSubsets|Extract a subset of sequences, given rules (skip X, print every X, max X, etc.)|
|excludeSeqsBySAM|Extract a subset of sequences, excluding those with no alignment in a SAM file|
|filterSeqSize|Extract a subset of sequences, retaining those within a range of lengths|
|getBiokangaAlignStats .pl|Produces alignment statistics from log-file output of 'biokanga align'|
|getSeqCGstats|Profiles sequences for GC% and base counts|
|getSeqCountTable|Produces a table counting occurances of individual sequences per inputs|
|getSeqQCStats|Profiles sequences for count, total/average/median bp length, av. GC%, fq phred scores...|
|getSeqSizeChart|Produces a table counting sequences of different lengths per inputs|
|getSeqSizeList|Print list of sequence bp lengths to stdout|
|getSeqSizeStats|Profiles sequences for total/average/median/min/max bp lengths|
|reverseComplement|Produces the reverse complements of sequences|
|splitInputs-snpTally-gz .pl|Companion to tallySNPs2, see README-tallySNPs.md|
|splitSeqsIntoXFiles|Will divide a sequence file up into multiple smaller sequence files|
|tallyGeneCoverageSamGZ|Produces a count of aligned reads per gene per sample|
|tallySNPs2|Counts aligned reads from different alleles at SNP positions, see README-tallySNPs.md|
|mergeKmerCounts|Merge Kmer count results from multiple samples into a multi-column table|

SeqReader.cpp/.h is a useful library for building upon.
It handles reading of fasta or fastq formatted sequence files and can handle .gz compressed inputs.
Allows for easy parsing with nextSeq() function and has various sequence manipulations built in. 

Code by Andrew Spriggs, CSIRO Ag&Food (www.csiro.au)
andrew.spriggs@csiro.au
github.com/spriggsy83

