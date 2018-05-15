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
#include "SNPTallyer2.h"
using namespace std;

/* By Andrew Spriggs, CSIRO Ag&Food, 2018 */
/* andrew.spriggs@csiro.au */
/* https://github.com/spriggsy83 */

/*See README-tallySNPs.txt*/

const char progName[] = "tallySNPs2";

bool getInputs(int argc, char* argv[], string& inRefSeqFileName, string& inSamplesFileName, string& outTabFilename, int& readDensityMin, int& edgeBuffer);
bool getSamples(const string& inSamplesFileName, vector<string>& labels, vector<string>& inSAMFileNames, vector<string>& inSNPFileNames);
void printHelp();

int main(int argc,char *argv[]){

	vector<string> inSAMFileNames;
	vector<string> inSNPFileNames;
	vector<string> labels;
	string inRefSeqFileName = "";
	string inSamplesFileName = "";
	string outTabFilename = "snpsOut.txt";
	int readDensityMin = 5;
	int edgeBuffer = 5;
	
	if(!getInputs(argc, argv, inRefSeqFileName, inSamplesFileName, outTabFilename, readDensityMin, edgeBuffer)){
		//cerr << "Process aborted.\n";
		return 1;
	}
	
	if(!getSamples(inSamplesFileName, labels, inSAMFileNames, inSNPFileNames)){
		cerr << "Process aborted.\n";
		return 1;
	}
	
	SNPTallyer theSNPTallyer(labels, inSAMFileNames, inSNPFileNames, inRefSeqFileName, outTabFilename, readDensityMin, edgeBuffer);
	
	if(theSNPTallyer.tallySNPs()){
		return 0;
	}else{
		return 1;
	}
}

bool getSamples(const string& inSamplesFileName, vector<string>& labels, vector<string>& inSAMFileNames, vector<string>& inSNPFileNames){
	
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
		if(lineParts.size() == 3){
			if(lineParts[0] != "" && lineParts[1] != "" && lineParts[2] != "" && 
					lineParts[0].find_first_of(' ') == string::npos && 
					lineParts[1].find_first_of(' ') == string::npos && 
					lineParts[2].find_first_of(' ') == string::npos){
				
				labels.push_back(lineParts[0]);		
				inSAMFileNames.push_back(lineParts[1]);
				inSNPFileNames.push_back(lineParts[2]);
			}
		}
	}
	infile.close();
	
	return true;
}

bool getInputs(int argc, char* argv[], string& inRefSeqFileName, string& inSamplesFileName, string& outTabFilename, int& readDensityMin, int& edgeBuffer){
	extern char *optarg;
	int opt;
	while ((opt = getopt(argc,argv,"i:r:o:d:e:h")) != EOF){
		switch(opt){
			case 'i':
				inSamplesFileName = optarg;
				break;
			case 'r':
				inRefSeqFileName = optarg;
				break;
			case 'o':
				outTabFilename = optarg;
				break;
			case 'd':
				readDensityMin = atoi(optarg);
				break;
			case 'e':
				edgeBuffer = atoi(optarg);
				break;
			case 'h':
			case '?':
			default:
				printHelp();
				return false;
		}
	}
	if(inSamplesFileName == "" || inRefSeqFileName == ""){
		printHelp();
		return false;
	}
	return true;
}

void printHelp(){
	cerr << "\t***** " << progName << " *****\n\t- Andrew Spriggs, CSIRO Ag&Food, 2018 -\n\n";
	cerr << "Usage:\t" << progName << " -i samplesFile -r refSeqFile [options]\n\n";
	cerr << "Options:\n";
	cerr << "\t-i samplesFile\t\tFilename of samples list (required) (see below)\n";
	cerr << "\t-r refSeqFile\t\tFilename of FASTA-formatted reference sequence (required)\n";
	cerr << "\t-o outTabFile\t\tFilename for table output (default = snpsOut.txt)\n";
	cerr << "\t-d readDepthMin\t\tMinimum read depth from a sample for a reported SNP (default = 5)\n";
	cerr << "\t-e edgeBuffer\t\tDon't count bases within __bp of ends of reads (default = 5)\n";
	cerr << "\n\n";
	cerr << "samplesFile should be a tab-separated file with each line representing a sample in the form:\n\n";
	cerr << "sample-name\tsam-file\tsnp-file\n\n";
	cerr << "...where sam-file is the filename for a SAM-formatted result of an alignment between the sample and the reference sequence, ";
	cerr << "snp-file is the Biokanga-Align-generated SNP-call csv file and sample-name is a short label to give the sample in outputs.\n\n";
	cerr << "Fasta, SAM and SNP files for input may be .gz compressed.\n\n";
}

