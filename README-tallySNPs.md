## TallySNPs
The tallySNPs program builds a table of counts of aligned reads from alternate alleles over a detected SNP, from multiple samples.
Input is assumed to be a series of 'biokanga align' results against the same reference, with SNP calling enabled.
(https://github.com/csiro-crop-informatics/biokanga)
Requires pairs of inputs from each sample- a SAM alignment file and a corresponding csv file of detected SNPs from Biokanga.
A SNP may be present for one sample but not others. Interrogation of SAM files builds details of exactly what reads/alleles are present in each.
E.g. answers 'is lack of SNP from lack of coverage or due to reads only matching reference?'

*tallySNP* re-reads each SAM input again for each reference sequence, to avoid holding all alignments in memory at once.
Thus it runs dramatically faster if inputs are split into smaller chunks beforehand.
This can be done with Perl script `splitInputs-snpTally-gz.pl` 

Example use-case:
```
#Create "alignList.txt", tab-separated, row per sample:
#sample-name  SAM-file.gz  snp-file.gz

REFSEQ=reference-sequences-aligned-to.fasta.gz
./getSeqSizeList $REFSEQ > refSizes.txt

mkdir -p splitInputs
mkdir -p splitOut

perl splitInputs-snpTally-gz.pl refSizes.txt $REFSEQ alignList.txt 500

for X in {0..499}
	do ./tallySNPs2 \
		-d 2 \
		-i splitInputs/fLists/inFilesList.pt${X}.txt \
		-r splitInputs/fasta/refSeqs.pt${X}.fasta.gz \
		-o splitOut/snpTally.pt${X}.txt
done

head -1 splitOut/snpTally.pt0.txt > snpTally2.result.txt
grep -v "RefID" splitOut/snpTally.pt* >> snpTally2.result.txt
perl -i -pe 's/^(\S+?)://' snpTally2.result.txt

rm -R splitInputs
rm -R splitOut
```