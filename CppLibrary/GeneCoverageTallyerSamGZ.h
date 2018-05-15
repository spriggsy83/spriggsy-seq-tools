#ifndef COVERTALLYER_H
#define COVERTALLYER_H

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
#include <algorithm>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
using namespace std;

/* By Andrew Spriggs, CSIRO Ag&Food, 2018 */
/* andrew.spriggs@csiro.au */
/* https://github.com/spriggsy83 */

struct GeneCoord {
	string name;
	unsigned int start;
	unsigned int end;

	GeneCoord(unsigned int newStart, unsigned int newEnd, const string& newName){
		name = newName;
		start = newStart;
		end = newEnd;
	}

	bool operator < (const GeneCoord& otherGene) const{
		return (start < otherGene.start);
	}
}; //!< gene start, gene end, gene name

typedef map< string, vector< GeneCoord > > CoordMap; //!< refSeqID, list of gene details
typedef map< string, vector< vector< unsigned int > > > ReadCountMap; //!< refSeqID, read counts per gene per sample

class GeneCoverageTallyer {
  private:
  	int minReads; //!< Minimum reads seen in any one sample to make it worth printing results for a coordinate range
	int updown; //!< Extra bases to include up/down-stream of gene coordinates
	vector<string> inSAMFileNames; //!< List of SAM alignment files for input
	vector<string> labels;  //!< List of sample names, one for each corresponding SAM input file
	string inCoordFileName; //!< Coordinate input file for tallying over
	int numSamples; //!< Number of input samples == number of SAM file inputs
	ofstream outtabfile; //!< Tabular output file name
	bool prepRan;  //!< Indicates that class has been initialised, output files have been opened and GeneCoverageTallyer is ready to run
	
	CoordMap geneCoords; //!< List of all starting coord positions, as read from input file
	ReadCountMap readCounts; //!< Read counts per refSeq, per gene, per sample
	unsigned int maxGene; //!< Largest gene size witnessed in gene list
	
		/*** Actual constructor code, called by constructor forms **/
	void prepareCoverageTallyer(const vector<string>& aLabelsList, const vector<string>& aSAMFileNamesList, const string& aCoordFileName, const string& aOutTabFileName, const int& aMinReads, const int& aUpDown);
		/*** Read coords list to form list of test coordinates.  Also initialises read count tables. **/
	bool loadCoordsList();
		/*** Read the sam.gz file to tally reads for a specific sample **/
	bool tallyReadsForSample(const int sNum);
		/*** Test if a line is SAM format aligned read then extract read coordinates **/
	bool getReadCoordFromSamLine(const string& line, unsigned int& rStart, unsigned int& rEnd, string& rRefID);
		/*** Test if a read overlaps a gene, add it to the tally **/
	void addReadTally(const unsigned int& rStart, const unsigned int& rEnd, const string& rRefID, const int sNum);
		/*** Finalise results to file **/
	bool writeOutput();
		
  public:
		/** Initialise with no defaults **/
	GeneCoverageTallyer(const vector<string>& aLabelsList, const vector<string>& aSAMFileNamesList, const string& aCoordFileName, const string& aOutTabFileName);
	GeneCoverageTallyer(const vector<string>& aLabelsList, const vector<string>& aSAMFileNamesList, const string& aCoordFileName, const string& aOutTabFileName, const int& aMinReads, const int& aUpDown);
	~GeneCoverageTallyer();
	
		/** Launch the full read tallying process **/
	bool tallyCoverage();
};

#endif

