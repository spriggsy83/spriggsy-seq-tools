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
#include "KmerCountMerger.h"
#include "SeqReader.h"
using namespace boost::iostreams;

/* By Andrew Spriggs, CSIRO Ag&Food, 2018 */
/* andrew.spriggs@csiro.au */
/* https://github.com/spriggsy83 */

/*** Initialise with defaults 
**/
KmerCountMerger::KmerCountMerger(const vector<string>& aLabelsList, const vector<string>& aFileNamesList, const string& aOutTabFileName){
	prepareKmerCountMerger(aLabelsList, aFileNamesList, aOutTabFileName, 20, false);
	return;
}

KmerCountMerger::KmerCountMerger(const vector<string>& aLabelsList, const vector<string>& aFileNamesList, const string& aOutTabFileName, const int& aMinCount, const bool& aTwoPassSet){
	prepareKmerCountMerger(aLabelsList, aFileNamesList, aOutTabFileName, aMinCount, aTwoPassSet);
	return;
}

/*** Actual constructor 
**/
void KmerCountMerger::prepareKmerCountMerger(const vector<string>& aLabelsList, const vector<string>& aFileNamesList, const string& aOutTabFileName, const int& aMinCount, const bool& aTwoPassSet){

	prepRan = false;
	labels = aLabelsList;
	inFileNames = aFileNamesList;
	minCount = aMinCount;
	twoPass = aTwoPassSet;
	numSamples = labels.size();
	
	outtabfile.open(aOutTabFileName.c_str());
	if(!outtabfile.is_open()){
		cerr << "Unable to open output file " << aOutTabFileName << "!\nNo processing will follow.\n";
		return;
	}

	prepRan = true;
	return;
}


KmerCountMerger::~KmerCountMerger(){
	if(outtabfile.is_open()){
		outtabfile.close();
	}
}

/*** Launch the full kmer count merging process
**/
bool KmerCountMerger::MergeKmerCounts(){
	
	if(!prepRan){
		return false;
	}
	
	for(int sNum=0; sNum < numSamples; sNum++){
		int fileType = testKmerFileForSample(sNum);
		switch(fileType){
			case 1:
				if(!readTabCountsForSample(sNum, 1)){
					cerr << "Failed to read kmers for " << labels[sNum] << " from tab file " << inFileNames[sNum] << endl;
				}
				break;
			case 2:
				if(!readFastaCountsForSample(sNum, 1)){
					cerr << "Failed to read kmers for " << labels[sNum] << " from fasta file" << inFileNames[sNum] << endl;
				}
				break;
			default:
				cerr << "Failed to read kmers for " << labels[sNum] << " from file " << inFileNames[sNum] << endl;
		}
	}
	if(twoPass){
		for(int sNum=0; sNum < numSamples; sNum++){
			int fileType = testKmerFileForSample(sNum);
			switch(fileType){
				case 1:
					if(!readTabCountsForSample(sNum, 2)){
						cerr << "Failed to read kmers for " << labels[sNum] << " from tab file " << inFileNames[sNum] << endl;
					}
					break;
				case 2:
					if(!readFastaCountsForSample(sNum, 2)){
						cerr << "Failed to read kmers for " << labels[sNum] << " from fasta file" << inFileNames[sNum] << endl;
					}
					break;
				default:
					cerr << "Failed to read kmers for " << labels[sNum] << " from file " << inFileNames[sNum] << endl;
			}
		}
	}

	if(!writeOutput()){
		return false;
	}
	return true;
}


/*** Peek at format of input file for a specific sample (0=error,1=tab,2=fasta)
**/
int KmerCountMerger::testKmerFileForSample(const int sNum){

	int fileType = 0;
	bool gzipFile = false;
	ifstream fileifs;
	string filename = inFileNames[sNum];
	if(filename.find("gz", filename.length()-3) != string::npos || 
			filename.find("GZ", filename.length()-3) != string::npos){
		gzipFile = true;
	}
	if(gzipFile){
		fileifs.open(filename.c_str(), ios_base::in | ios_base::binary);
	}else{
		fileifs.open(filename.c_str(), ios_base::in );
	}
	if (!fileifs.is_open()){
		cerr << "Unable to open file " << filename << "!\n";
		return false;
	}else if (!fileifs.good()){
		cerr << "File " << filename << " is empty!\n";
		fileifs.close();
		return false;
	}
	try {
		filtering_istream infile;
		if(gzipFile){
			infile.push(gzip_decompressor());
		}
		infile.push(fileifs);
		string line;
		getline(infile, line);
		if(line.find('>') == 0){
			fileType = 2;
		}else if(line.find('\t') > 0 && line.find('\t') != std::string::npos){
			fileType = 1;
		}
		fileifs.close();
	}
	catch(const gzip_error& e) {
		cerr << "Error while reading file " << filename << endl;
		cerr << e.what() << endl;
		fileifs.close();
		return 0;
	}
	return fileType;
}


/*** Read the input file for a specific sample - tab format 
**/
bool KmerCountMerger::readTabCountsForSample(const int sNum, const int pass){
	ifstream fileifs;
	bool gzipFile = false;
	string filename = inFileNames[sNum];
	if(filename.find("gz", filename.length()-3) != string::npos || 
			filename.find("GZ", filename.length()-3) != string::npos){
		gzipFile = true;
	}
	if(gzipFile){
		fileifs.open(filename.c_str(), ios_base::in | ios_base::binary);
	}else{
		fileifs.open(filename.c_str(), ios_base::in );
	}
	if (!fileifs.is_open()){
		cerr << "Unable to open file " << filename << "!\n";
		return false;
	}else if (!fileifs.good()){
		cerr << "File " << filename << " is empty!\n";
		fileifs.close();
		return false;
	}
	try {
		filtering_istream infile;
		if(gzipFile){
			infile.push(gzip_decompressor());
		}
		infile.push(fileifs);
		cout << "Parsing kmer counts from " << filename << endl;

		unsigned int totKmers = 0;

		string line;
		while(getline(infile, line)){
			stringstream linestream(line);
			vector<string> lineParts;
			lineParts.reserve(2);
			string aLinePart;
			// Tab separated split
			while(getline(linestream, aLinePart, '\t')){
				lineParts.push_back(aLinePart);
			}
			if(lineParts.size() == 2){
				string kmerSeq = lineParts[0];
				unsigned int kmerCount;
				stringstream valuess(lineParts[1]);
				valuess >> kmerCount;

				if( (twoPass && pass == 1 && kmerCount >= minCount) || !twoPass ){
					totKmers++;
					if(kmerCounts.count(kmerSeq) == 0){
						kmerCounts[kmerSeq] = vector< unsigned int >(numSamples, 0);
					}
					kmerCounts[kmerSeq][sNum] = kmerCount;
				}else if(twoPass && pass == 2 && kmerCount < minCount){
					// Second pass, add sample count if found with minCount in another sample
					if(kmerCounts.count(kmerSeq) == 1){
						totKmers++;
						kmerCounts[kmerSeq][sNum] = kmerCount;
					}
				}
			}
		}
		fileifs.close();
		cout << "Loaded " << totKmers << " kmer seqs and counts from " << filename << endl;
	}
	catch(const gzip_error& e) {
		cerr << "Error while reading file " << filename << endl;
		cerr << e.what() << endl;
		fileifs.close();
		return false;
	}
	return true;
}


/*** Read the input file for a specific sample - FASTA format 
**/
bool KmerCountMerger::readFastaCountsForSample(const int sNum, const int pass){
	string filename = inFileNames[sNum];
	try {
		SeqReader inFile(filename);
		unsigned int totKmers = 0;
		while(inFile.nextSeq()){
			unsigned int kmerCount;
			stringstream valuess(inFile.getSeqID());
			valuess >> kmerCount;
			if(inFile.getSeqLen() > 0){
				totKmers++;
				string kmerSeq = inFile.getSeq();
				if( (twoPass && pass == 1 && kmerCount >= minCount) || !twoPass ){
					totKmers++;
					if(kmerCounts.count(kmerSeq) == 0){
						kmerCounts[kmerSeq] = vector< unsigned int >(numSamples, 0);
					}
					kmerCounts[kmerSeq][sNum] = kmerCount;
				}else if(twoPass && pass == 2 && kmerCount < minCount){
					// Second pass, add sample count if found with minCount in another sample
					if(kmerCounts.count(kmerSeq) == 1){
						totKmers++;
						kmerCounts[kmerSeq][sNum] = kmerCount;
					}
				}
			}
		}
		cout << "Loaded " << totKmers << " kmer seqs and counts from " << filename << endl;
	}
	catch(const gzip_error& e) {
		cerr << "Error while reading file " << filename << endl;
		cerr << e.what() << endl;
		return false;
	}
	return true;
}


/*** Finalise results to file
**/
bool KmerCountMerger::writeOutput(){

	if(!outtabfile.is_open()){
		cerr << "Output file is not open for writing! Can't output results.\n";
		return false;
	}
	
	outtabfile << "KmerSeq";
	for(int sNum=0; sNum < numSamples; sNum++){
		outtabfile << "\t" << labels[sNum];
	}
	outtabfile << "\n";

	unsigned int totKmers = 0;
	for(KmerCountMap::iterator aKmer=kmerCounts.begin(); aKmer!=kmerCounts.end(); ++aKmer){
		string kmerSeq = aKmer->first;
		vector< unsigned int > countsV = aKmer->second;

		// Check if minimum count met for printing
		bool doPrint = false;
		for(int sNum=0; sNum < numSamples && !doPrint; sNum++){
			if(countsV[sNum] >= minCount){
				doPrint = true;
			}
		}

		// Output for a kmer
		if(doPrint){
			totKmers++;
			outtabfile << kmerSeq;
			for(int sNum=0; sNum < numSamples; sNum++){
				outtabfile << "\t" << countsV[sNum];
			}
			outtabfile << "\n";
		}
	}
	outtabfile.close();
	cout << "Wrote " << totKmers << " kmer seqs and counts to out file" << endl;
}

