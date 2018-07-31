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
#include "GeneCoverageTallyerSamGZ.h"
using namespace std;

/* By Andrew Spriggs, CSIRO Ag&Food, 2018 */
/* andrew.spriggs@csiro.au */
/* https://github.com/spriggsy83 */

/*** Initialise with defaults 
**/
GeneCoverageTallyer::GeneCoverageTallyer(const vector<string>& aLabelsList, const vector<string>& aSAMFileNamesList, const string& aCoordFileName, const string& aOutTabFileName){
	prepareCoverageTallyer(aLabelsList, aSAMFileNamesList, aCoordFileName, aOutTabFileName, 20, 0);
	return;
}

GeneCoverageTallyer::GeneCoverageTallyer(const vector<string>& aLabelsList, const vector<string>& aSAMFileNamesList, const string& aCoordFileName, const string& aOutTabFileName, const int& aMinReads, const int& aUpDown){
	prepareCoverageTallyer(aLabelsList, aSAMFileNamesList, aCoordFileName, aOutTabFileName, aMinReads, aUpDown);
	return;
}

/*** Actual constructor 
**/
void GeneCoverageTallyer::prepareCoverageTallyer(const vector<string>& aLabelsList, const vector<string>& aSAMFileNamesList, const string& aCoordFileName, const string& aOutTabFileName, const int& aMinReads, const int& aUpDown){

	prepRan = false;
	labels = aLabelsList;
	inSAMFileNames = aSAMFileNamesList;
	inCoordFileName = aCoordFileName;
	minReads = aMinReads;
	updown = aUpDown;
	numSamples = labels.size();
	
	outtabfile.open(aOutTabFileName.c_str());
	if(!outtabfile.is_open()){
		cerr << "Unable to open output file " << aOutTabFileName << "!\nNo read tallying will follow.\n";
		return;
	}

	prepRan = true;
	return;
}


GeneCoverageTallyer::~GeneCoverageTallyer(){
	if(outtabfile.is_open()){
		outtabfile.close();
	}
}

/*** Launch the full read tallying process
**/
bool GeneCoverageTallyer::tallyCoverage(){
	
	if(!prepRan){
		return false;
	}
	
	if(!loadCoordsList()){
		cerr << "Failed to load any input test coordinates." << endl;
		return false;
	}else{
		
		#pragma omp parallel for
		for(int sNum=0; sNum < numSamples; sNum++){	
			if(!tallyReadsForSample(sNum)){
				cerr << "Failed to parse reads for " << labels[sNum] << " from " << inSAMFileNames[sNum] << endl;
			}
		}

		if(!writeOutput()){
			return false;
		}
	}
	return true;
}


/*** Read gene coords list to form list of test coordinates
** Also initialises read count tables
**/
bool GeneCoverageTallyer::loadCoordsList(){
	ifstream infile;
	string line;
	infile.open(inCoordFileName.c_str());
	if (!infile.is_open()){
		cerr << "Unable to open coords file " << inCoordFileName << "!\n";
		return false;
	}else if (!infile.good()){
		cerr << "Coords file " << inCoordFileName << " is empty!\n";
		infile.close();
		return false;
	}
	
	cout << "Parsing coords list from " << inCoordFileName << endl;

	maxGene = 0;
	
	while(getline(infile, line)){
		stringstream linestream(line);
		vector<string> lineParts;
		lineParts.reserve(4);
		string aLinePart;
		// Tab separated split
		while(getline(linestream, aLinePart, '\t')){
			lineParts.push_back(aLinePart);
		}
		
		//refID \t start \t end \t more
		if(lineParts.size() >= 4){
			string refID = lineParts[0];
			string geneName = lineParts[3];
			stringstream startSS(lineParts[1]);
			unsigned int start;
			startSS >> start;
			stringstream endSS(lineParts[2]);
			unsigned int end;
			endSS >> end;
			if(end < start){
				unsigned int temp = start;
				start = end;
				end = temp;
			}

			if(geneCoords.count(refID) == 0){
				geneCoords[refID] = vector< GeneCoord >();
				readCounts[refID] = vector< vector< unsigned int > >();
			}
			
			geneCoords[refID].push_back( GeneCoord(start, end, geneName) );
			readCounts[refID].push_back( vector< unsigned int >(numSamples, 0) );

			if(end - start + 1 > maxGene){
				maxGene = end - start + 1;
			}
		}
	}
	infile.close();
	
	int numRefIDs = geneCoords.size();
	unsigned int numCoords = 0;
	for (CoordMap::iterator aRef=geneCoords.begin(); aRef!=geneCoords.end(); ++aRef){
		sort(aRef->second.begin(), aRef->second.end());
		numCoords += aRef->second.size();
	}
	cout << "Loaded " << numCoords << " test coords over " << numRefIDs << " reference sequences." << endl;
	
	if(numCoords == 0){
		return false;
	}
	return true;
}

/*** Read the sam.gz file to tally reads for a specific sample
**/
bool GeneCoverageTallyer::tallyReadsForSample(const int sNum){

	ifstream fileifs(inSAMFileNames[sNum].c_str(), ios_base::in | ios_base::binary);
	try {
		boost::iostreams::filtering_istream infile;
		infile.push(boost::iostreams::gzip_decompressor());
		infile.push(fileifs);
		cout << "Parsing reads from " << inSAMFileNames[sNum] << endl;
		unsigned int sampTotReads = 0;
		string line;
		while(getline(infile, line)){
			unsigned int rStart = 0;
			unsigned int rEnd = 0;
			string rRefID;
			if(getReadCoordFromSamLine(line, rStart, rEnd, rRefID)){
				sampTotReads++;
				addReadTally(rStart, rEnd, rRefID, sNum);
			}
		}
		cout << "Loaded " << sampTotReads << " reads from " << inSAMFileNames[sNum] << endl;
	}
	catch(const boost::iostreams::gzip_error& e) {
		cerr << "Error while reading sam.gz file " << inSAMFileNames[sNum] << endl;
		cerr << e.what() << endl;
		return false;
	}
	return true;
}


/*** Test if a line is SAM format aligned read then extract read coordinates
**/
bool GeneCoverageTallyer::getReadCoordFromSamLine(const string& line, unsigned int& rStart, unsigned int& rEnd, string& rRefID ){

	stringstream linestream(line);
	vector<string> lineParts;
	lineParts.reserve(11);
	string aLinePart;
	// Tab separated split
	while(getline(linestream, aLinePart, '\t')){
		lineParts.push_back(aLinePart);
	}
	if(lineParts.size() == 11){
		// readID == [0], refID == [2], start == [3], cigar == [5], readSeq == [9]

		// If there are genes to tally for refSeq aligned to
		if(geneCoords.count(lineParts[2])){
			rRefID = lineParts[2];

			stringstream startstream(lineParts[3]);
			startstream >> rStart;
			rStart = rStart - 1;
			rEnd = rStart;

			// Process cigar string
			int valStartI = 0;
			for(int i=0; i<lineParts[5].length(); i++){
				if(!isdigit(lineParts[5].at(i))){
					
					char action = lineParts[5].at(i);
					string valuestr = lineParts[5].substr(valStartI, i-valStartI);
					stringstream valuess(valuestr);
					int value;
					valuess >> value;
					
					switch (action){
						case('M'):  // Match
						case('='):  // Match
						case('X'):  // Match
							rEnd += value;
							break;
							
						case('S'):  // Soft-trim
							break;
							
						case('N'): // Intron
							rEnd += value;
							break;
							
						case('D'):  // Del
							rEnd += value;
							break;
							
						case('I'):  // Ins
							break;
					}
					valStartI = i+1;
				}
			}

			if(rEnd > rStart){
				return true;
			}
		}else{
			return false;
		}
	}else{
		return false;
	}
}


/*** Test if a read overlaps a gene or genes, add it to the tally 
**/
void GeneCoverageTallyer::addReadTally(const unsigned int& rStart, const unsigned int& rEnd, const string& rRefID, const int sNum){

	// Binary search to first potentially overlapping gene
	int lowBound = 0;
	int highBound = geneCoords[rRefID].size();
	while(lowBound != highBound){
		int midpoint = (lowBound + highBound) / 2;
		if( geneCoords[rRefID][midpoint].start + maxGene + updown <= rStart ){
			lowBound = midpoint + 1;
		}else{
			highBound = midpoint;
		}
	}
	for(int i=lowBound; i < geneCoords[rRefID].size(); i++ ){

		// If gene and read overlap, including up/down-stream buffer
		if( (rStart + updown >= geneCoords[rRefID][i].start && rStart < geneCoords[rRefID][i].end + updown) || 
				(rEnd + updown > geneCoords[rRefID][i].start && rEnd <= geneCoords[rRefID][i].end + updown) || 
				(rStart + updown < geneCoords[rRefID][i].start && rEnd > geneCoords[rRefID][i].end + updown) ){

			// Add read to gene tally for sample
			readCounts[rRefID][i][sNum] += 1;
		}

		if(geneCoords[rRefID][i].start > rEnd + updown){
			break;
		}
	}
	return;
}


/*** Finalise results to file
**/
bool GeneCoverageTallyer::writeOutput(){

	if(!outtabfile.is_open()){
		cerr << "Output file is not open for writing! Can't output results.\n";
		return false;
	}
	
	outtabfile << "Gene\tRefID\tStart\tEnd";
	for(int sNum=0; sNum < numSamples; sNum++){
		outtabfile << "\t" << labels[sNum];
	}
	outtabfile << "\n";

	for(CoordMap::iterator aRef=geneCoords.begin(); aRef!=geneCoords.end(); ++aRef){
		string refID = aRef->first;

		for(int geneI = 0; geneI < geneCoords[refID].size(); geneI++ ){

			// Check if minimum reads met for printing
			bool doPrint = false;
			for(int sNum=0; sNum < numSamples && !doPrint; sNum++){
				if(readCounts[refID][geneI][sNum] >= minReads){
					doPrint = true;
				}
			}

			// Output for a gene
			if(doPrint){
				outtabfile << geneCoords[refID][geneI].name << "\t";
				outtabfile << refID << "\t";
				outtabfile << geneCoords[refID][geneI].start << "\t";
				outtabfile << geneCoords[refID][geneI].end;
				for(int sNum=0; sNum < numSamples; sNum++){
					outtabfile << "\t" << readCounts[refID][geneI][sNum];
				}
				outtabfile << "\n";
			}
		}
	}

	outtabfile.close();
}

