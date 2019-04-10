#ifndef KMERTALLYER_H
#define KMERTALLYER_H

#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <string>
#include <cstring>
#include <ctype.h>
#include <sstream>
#include <algorithm>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include "SeqReader.h"
using namespace std;

/* By Andrew Spriggs, CSIRO Ag&Food, 2018 */
/* andrew.spriggs@csiro.au */
/* https://github.com/spriggsy83 */

typedef map< string, vector< unsigned int > > KmerCountMap; //!< kmer seq, read counts per kmer seq per sample

class KmerCountMerger {
  private:
  	int minCount; //!< Minimum reads seen in any one sample to make it worth printing results for a kmer
	bool twoPass; //!< Two-pass mode on/off
	vector<string> inFileNames; //!< List of kmer-count tab-sep/fasta files for input
	vector<string> labels;  //!< List of sample names, one for each corresponding input file
	int numSamples; //!< Number of input samples == number of tab file inputs
	ofstream outtabfile; //!< Tabular output file name
	bool prepRan;  //!< Indicates that class has been initialised, output files have been opened and processing is ready to run
	
	KmerCountMap kmerCounts; //!< Read counts per kmer seq per sample
	
		/*** Actual constructor code, called by constructor forms **/
	void prepareKmerCountMerger(const vector<string>& aLabelsList, const vector<string>& aFileNamesList, const string& aOutTabFileName, const int& aMinCount, const bool& aTwoPassSet);
		/*** Peek at format of input file for a specific sample (0=error,1=tab,2=fasta) **/
	int testKmerFileForSample(const int sNum);
		/*** Read the input file for a specific sample - tab format **/
	bool readTabCountsForSample(const int sNum, const int pass);
		/*** Read the input file for a specific sample - FASTA format **/
	bool readFastaCountsForSample(const int sNum, const int pass);
		/*** Finalise results to file **/
	bool writeOutput();

  public:
		/** Initialise with no defaults **/
	KmerCountMerger(const vector<string>& aLabelsList, const vector<string>& aFileNamesList, const string& aOutTabFileName);
	KmerCountMerger(const vector<string>& aLabelsList, const vector<string>& aFileNamesList, const string& aOutTabFileName, const int& aMinCount, const bool& aTwoPassSet);
	~KmerCountMerger();
	
		/** Launch the full kmer count merging process **/
	bool MergeKmerCounts();
};

#endif
