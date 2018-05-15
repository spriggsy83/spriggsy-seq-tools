#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include "SeqReader.h"
using namespace std;

/* By Andrew Spriggs, CSIRO Ag&Food, 2018 */
/* andrew.spriggs@csiro.au */
/* https://github.com/spriggsy83 */

const char progName[] = "filterSeqSize";

bool getInputs(int argc, char* argv[], string& inFileName, string& outFileName, int& minSize, int& maxSize);

int main(int argc,char *argv[]){

	string inFileName;
	string outFileName;
	int minSize = -1;
	int maxSize = -1;
	int countKept = 0;
	int countReject = 0;
	
	if(!getInputs(argc, argv, inFileName, outFileName, minSize, maxSize)){
		cerr << "Process aborted.\n";
		return 0;
	}
	
	ofstream outfile(outFileName.c_str());
	if(!outfile.is_open()){
		cerr << "Unable to open output file " << outFileName << "!\nProcess aborted.\n";
		return 0;
	}
	
	SeqReader inFile(inFileName);
	
	while(inFile.nextSeq()){
		int len = inFile.getSeqLen();
		
		if((minSize == -1 || len >= minSize) && (maxSize == -1 || len <= maxSize)){
			outfile << inFile.toString();
			countKept++;
		}else{
			countReject++;
		}
	}
	
	cout << "Filtered " << inFileName << " for min:" << minSize << " to max:" << maxSize << "\n";
	cout << "Retained = " << countKept << " | Rejected = " << countReject << "\n";
	return 0;
}

bool getInputs(int argc, char* argv[], string& inFileName, string& outFileName, int& minSize, int& maxSize){
	if(argc != 5){
		cerr << "\t***** " << progName << " *****\n\t- Andrew Spriggs, CSIRO Ag&Food, 2018 -\n";
		cerr << "Extract subset of sequences within a length range.\n";
		cerr << "Input file may be fasta or fastq and may be .gz compressed.\n";
		cerr << "Correct command line usage:\n" << argv[0] << " <infile> <outfile> <min seq len> <max seq len>\n";
		cerr << "Give a min or max of -1 for no limit\nMax and min are inclusive\n";
		return false;
	}
	
	inFileName = argv[1];
	outFileName = argv[2];
	minSize = atoi(argv[3]);
	maxSize = atoi(argv[4]);
	
	return true;
}

