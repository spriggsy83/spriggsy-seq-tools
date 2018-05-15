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
#include "GeneCoverageTallyerSamGZ.h"
using namespace std;

/* By Andrew Spriggs, CSIRO Ag&Food, 2018 */
/* andrew.spriggs@csiro.au */
/* https://github.com/spriggsy83 */

const char progName[] = "tallyGeneCoverageSamGZ";

bool getInputs(int argc, char* argv[], string& inSamplesFileName, string& inCoordsFileName, string& outTabFilename, int& minReads, int& updown);
bool getSamples(const string& inSamplesFileName, vector<string>& labels, vector<string>& inSAMFileNames);
void printHelp();

int main(int argc,char *argv[]){

	vector<string> inSAMFileNames;
	vector<string> labels;
	string inCoordFileName = "";
	string inSamplesFileName = "";
	string outTabFilename = "";
	int minReads = 20;
	int updown = 0;
	
	if(!getInputs(argc, argv, inSamplesFileName, inCoordFileName, outTabFilename, minReads, updown)){
		//cerr << "Process aborted.\n";
		return 1;
	}
	
	if(!getSamples(inSamplesFileName, labels, inSAMFileNames)){
		cerr << "Process aborted.\n";
		return 1;
	}
	
	GeneCoverageTallyer theCoverageTallyer(labels, inSAMFileNames, inCoordFileName, outTabFilename, minReads, updown);
	
	if(theCoverageTallyer.tallyCoverage()){
		return 0;
	}else{
		return 1;
	}
}

bool getSamples(const string& inSamplesFileName, vector<string>& labels, vector<string>& inSAMFileNames){
	
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
				inSAMFileNames.push_back(lineParts[1]);
			}
		}
	}
	infile.close();
	
	return true;
}

bool getInputs(int argc, char* argv[], string& inSamplesFileName, string& inCoordsFileName, string& outTabFilename, int& minReads, int& updown){
	extern char *optarg;
	int opt;
	while ((opt = getopt(argc,argv,"i:c:o:m:u:h")) != EOF){
		switch(opt){
			case 'i':
				inSamplesFileName = optarg;
				break;
			case 'c':
				inCoordsFileName = optarg;
				break;
			case 'o':
				outTabFilename = optarg;
				break;
			case 'm':
				minReads = atoi( optarg );
				break;
			case 'u':
				updown = atoi( optarg );
				break;
			case 'h':
			case '?':
			default:
				printHelp();
				return false;
		}
	}
	if(inSamplesFileName == "" || inCoordsFileName == "" || outTabFilename == ""){
		printHelp();
		return false;
	}
	return true;
}

void printHelp(){
	cerr << "\t***** " << progName << " *****\n\t- Andrew Spriggs, CSIRO Ag&Food, 2018 -\n";
	cerr << "Produces a count of aligned reads per gene per sample.\n\n";
	cerr << "Usage:\t" << progName << " -i samplesFile -c coordsFile -o outTabFile\n\n";
	cerr << "Parameters:\n";
	cerr << "\t-i samplesFile\t\tFilename of samples list\n";
	cerr << "\t-c coordsFile\t\tFilename of gene coordinates to test coverage over\n";
	cerr << "\t-o outTabFile\t\tFilename for read tally table output\n";
	cerr << "\t-m minReads\t\tMinimum reads from a sample required for results over a coordinate range to be printed (default=20)\n";
	cerr << "\t-u updown\t\tExtra bases to include up/down-stream of gene coordinates (default=0)\n";
	cerr << "\n\n";
	cerr << "samplesFile should be a tab-separated file with each line representing a sample in the form:\n\n";
	cerr << "sample-name\tsam.gz-file\n\n";
	cerr << "...where sam.gz-file is the filename for a GZipped SAM-formatted result of an alignment between the sample and the reference sequence";
	cerr << "and sample-name is a short label to give the sample in outputs.\n";
	cerr << "...and coordsFile is the filename of a tab-separated file containing coords to test coverage over, in form-\n";
	cerr << "refID\tstart-coord\tend-coord\tgene-name\t[additional columns]\n\n";
}

