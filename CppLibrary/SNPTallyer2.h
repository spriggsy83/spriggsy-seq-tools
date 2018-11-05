#ifndef SNPTALLYER_H
#define SNPTALLYER_H

#include <omp.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <string>
#include <cstring>
#include <ctype.h>
#include <sstream>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include "SeqReader.h"
#include "AlignedRead.h"
using namespace std;

/* By Andrew Spriggs, CSIRO Ag&Food, 2018 */
/* andrew.spriggs@csiro.au */
/* https://github.com/spriggsy83 */

/*** Prediction of SNPs from sequencing of tetraploid genomes, using a single sub-genome as reference
**/
class SNPTallyer {
  private:
		// Initialisations for defaultRefLabel and defaultAltRefLabel in SNPTallyer.cpp
	static const int defaultReadDepthMin = 5; //!< Default readDepthMin =5
	static const int defaultEdgeBuffer = 5; //!< Default edgeBuffer =5
	static const int minorAlleleThresh = 4; //!< SNP looks real if a minor allele has less than 1/n reads of SNP allele
	static const int readEndBuffer = 500; //!< Offset reads by this when finding search start in tallyBases algorithm
	
	vector<string> inSAMFileNames; //!< List of SAM alignment files for input
	vector<string> inSNPFileNames; //!< List of biokanga-align SNP files for input
	vector<string> labels;  //!< List of sample names, one for each corresponding SAM input file
	string inRefSeqFileName;  //!< Fasta reference sequence file name for input
	int numSamples; //!< Number of input samples == number of SAM file inputs
	vector<omp_lock_t> ompWriteLocks; //!< Write-locks for parellelisation. Multi for loops-within-loops
	ofstream outtabfile;  //!< Tabular output file name
	int outFormat;  //!< Output format, 1 = row per SNP, 2 = row per allele, 3 = row per pos
	bool filesReady;  //!< Indicates that class has been initialised, output files have been opened and SNPTallyer is ready to run
	bool snpsPreLoaded;  //!< Indicates that biokanga-align SNP files have been parsed and snpPreList prepared
	string readsLoadedFor;	//!< Stores reference sequence ID of currently loaded read alignments
	int readDepthMin; //!< Minimum read depth from a sample for a reported SNP 
	int edgeBuffer; //!< In test of reads spanning SNPs, this adds an untested buffer to edge of read
	map< string, set<unsigned int> > snpPreList; //!< List of all starting SNP positions, as read from the biokanga-align SNP files
	vector<vector<AlignedRead> > reads; //!< Stores details of aligned reads for a reference sequence
	
		/*** Actual constructor code, called by constructor forms **/
	void prepareSNPTallyer(const vector<string>& aLabelsList, 
						const vector<string>& aSAMFileNamesList, 
						const vector<string>& aSNPFileNamesList, 
						const string& aRefSeqFileName, 
						const string& aOutTabFileName, 
						const int aOutFormat, 
						const int aReadDepthMin, 
						const int aEdgeBuffer);

		/*** Read Biokanga-Align SNP lists to form starting list of SNP locations **/
	bool loadSNPLists();

		/*** Load read alignments against a reference seq, from SAM files, for all samples **/
	void readReadsAll(const string& refID, int refSeqLen);

		/*** Load read alignments against a reference seq, from SAM files, for a single sample, adding to vector of aligned reads **/
	int readReadsSample(const string& inSamFileName, 
						string refID, 
						vector<AlignedRead>& reads);

		/*** Test reads from all samples over a SNP coord and print results if suitable **/
	bool testSNP(const unsigned int snpCoord, 
				const char refBase, 
				const string& refID);

		/*** Tally bases from reads seen over a SNP coord for a sample **/
	unsigned int tallyBases(const unsigned int snpCoord, 
							const vector<AlignedRead>& sampReads, 
							vector<unsigned int>& baseTally);
		
  public:
		/** Initialise with no defaults **/
	SNPTallyer(const vector<string>& aLabelsList, 
				const vector<string>& aSAMFileNamesList, 
				const vector<string>& aSNPFileNamesList, 
				const string& aRefSeqFileName, 
				const string& aOutTabFileName, 
				const int aOutFormat, 
				const int aReadDepthMin, 
				const int aEdgeBuffer);
	~SNPTallyer();
	
		/** Launch SNP tally across all reference sequences **/
	bool tallySNPs();
		/** Launch SNP detection against a single reference sequence **/ 
	bool tallySNPsOnRef(const string& refSeq, const string& refID);
};

#endif

