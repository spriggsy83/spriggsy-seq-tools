#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <unistd.h>
#include "KmerCountMerger.h"
using namespace std;

/* By Andrew Spriggs, CSIRO Ag&Food, 2018 */
/* andrew.spriggs@csiro.au */
/* https://github.com/spriggsy83 */

bool getInputs(int argc, char* argv[], string& inSamplesFileName, string& outTabFilename, int& minCount, bool& twoPass);
bool getSamples(const string& inSamplesFileName, vector<string>& labels, vector<string>& inFileNames);
void printHelp();

int main(int argc,char *argv[]){

	vector<string> inFileNames;
	vector<string> labels;
	string inSamplesFileName = "";
	string outTabFilename = "";
	int minCount = 20;
	bool twoPass = false;
	
	if(!getInputs(argc, argv, inSamplesFileName, outTabFilename, minCount, twoPass)){
		//cerr << "Process aborted.\n";
		return 1;
	}
	
	if(!getSamples(inSamplesFileName, labels, inFileNames)){
		cerr << "Process aborted.\n";
		return 1;
	}
	
	KmerCountMerger theKmerCountMerger(labels, inFileNames, outTabFilename, minCount, twoPass);
	
	if(theKmerCountMerger.MergeKmerCounts()){
		return 0;
	}else{
		return 1;
	}
}

bool getSamples(const string& inSamplesFileName, vector<string>& labels, vector<string>& inFileNames){
	
	ifstream infile;
	string line;
	infile.open(inSamplesFileName.c_str());
	if (!infile.is_open()){
		cerr << "Unable to open samples file " << inSamplesFileName << "!\n";
		return false;
	}else if (!infile.good()){
		cerr << "Samples file " << inSamplesFileName << " is empty!\n";
		infile.close();
		return false;
	}
	
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
			if(lineParts[0] != "" && lineParts[1] != "" && 
					lineParts[0].find_first_of(' ') == string::npos && 
					lineParts[1].find_first_of(' ') == string::npos){
				
				labels.push_back(lineParts[0]);		
				inFileNames.push_back(lineParts[1]);
			}
		}
	}
	infile.close();
	
	return true;
}

bool getInputs(int argc, char* argv[], string& inSamplesFileName, string& outTabFilename, int& minCount, bool& twoPass){
	twoPass = false;
	extern char *optarg;
	int opt;
	while ((opt = getopt(argc,argv,"i:o:m:2h")) != EOF){
		switch(opt){
			case 'i':
				inSamplesFileName = optarg;
				break;
			case 'o':
				outTabFilename = optarg;
				break;
			case 'm':
				minCount = atoi( optarg );
				break;
			case '2':
				twoPass = true;
				break;
			case 'h':
			case '?':
			default:
				printHelp();
				return false;
		}
	}
	if(inSamplesFileName == "" || outTabFilename == ""){
		printHelp();
		return false;
	}
	return true;
}

void printHelp(){
	cerr << "\t***** mergeKmerCountsFrom *****\n\t- Andrew Spriggs, CSIRO Ag&Food, 2018 -\n\n";
	cerr << "Usage:\tmergeKmerCountsFrom -i samplesFile -o outTabFile\n\n";
	cerr << "Parameters:\n";
	cerr << "\t-i samplesFile\t\tFilename of samples list\n";
	cerr << "\t-o outTabFile\t\tFilename for kmer counts table output\n";
	cerr << "\t-m minCount\t\tMinimum count of a kmer from any sample required for kmer to be printed (default=20)\n";
	cerr << "\t-2\t\tEnable two-pass mode, which may be faster with higher minCounts\n";
	cerr << "\n\n";
	cerr << "samplesFile should be a tab-separated file with each line representing a sample in the form:\n";
	cerr << "sample-name\tkmer-count-file\n\n";
	cerr << "...where kmer-count-file is the filename for either\n";
	cerr << "1. a tab-separated result of kmer counting (kmer[tab]count)\n";
	cerr << "OR 2. a FASTA-formatted result of kmer counting (seq ID as count),\n";
	cerr << "and sample-name is a short label to give the sample in outputs.\n";
	cerr << "Kmer count files may be .gz compressed.\n\n";
}
