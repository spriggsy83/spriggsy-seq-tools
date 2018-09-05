#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include "SeqReader.h"
using namespace std;

/* By Andrew Spriggs, CSIRO Ag&Food, 2018 */
/* andrew.spriggs@csiro.au */
/* https://github.com/spriggsy83 */

const char progName[] = "extractSeqSubsets";

void printHelp();
bool getInputs(int argc, char* argv[], string& inFileName, string& outFileName, 
			int& skipNum, unsigned long& maxNum, int& maxGbp, int& mode, int& X);

int main(int argc,char *argv[]){

	string inFileName;
	string outFileName;
	int skipNum = 0;
	unsigned long maxNum = 0;
	int maxGbp = 0;
	unsigned long maxbp = 0;
	int mode = 0; // 0 = print-all, 1 = extract-every-X mode, 2 = exclude-every-X mode
	int X = 0;
	
	if(!getInputs(argc, argv, inFileName, outFileName, skipNum, maxNum, maxGbp, mode, X)){
		cerr << "Process aborted.\n";
		return 0;
	}
	maxbp = (unsigned long)maxGbp * 1000000000;
	
	ofstream outfile(outFileName.c_str());
	if(!outfile.is_open()){
		cerr << "Unable to open output file " << outFileName << "!\nProcess aborted.\n";
		return 0;
	}
	
	SeqReader inFile(inFileName);
	
	int counter = 0;
	int skipCount = 0;
	unsigned long printbp = 0;
	unsigned long numPrinted = 0;

	while( inFile.nextSeq() && (numPrinted < maxNum || maxNum == 0) && (printbp < maxbp || maxbp == 0) ){

		if(skipNum > 0 && skipCount < skipNum){
			skipCount++;

		}else if(mode == 0){  // Print all
				outfile << inFile.toString();
				printbp += inFile.getSeqLen();
				numPrinted++;

		}else if(mode == 1){  // extract-every-X mode
			counter++;
			if(counter == X){
				counter = 0;
				outfile << inFile.toString();
				printbp += inFile.getSeqLen();
				numPrinted++;
			}
			
		}else if(mode == 2){  // exclude-every-X mode
			counter++;
			if(counter == X){
				counter = 0;
			}else{
				outfile << inFile.toString();
				printbp += inFile.getSeqLen();
				numPrinted++;
			}
		}
	}
	
	cout << "Output " << numPrinted << "sequences (" << printbp << " bp).\n";
	return 0;
}

bool getInputs(int argc, char* argv[], string& inFileName, string& outFileName, 
			int& skipNum, unsigned long& maxNum, int& maxGbp, int& mode, int& X){
	extern char *optarg;
	int opt;
	mode = 0; // 0 = print-all, 1 = extract-every-X mode, 2 = exclude-every-X mode
	X = 0;
	skipNum = 0;
	maxNum = 0;
	maxGbp = 0;
	while ((opt = getopt(argc,argv,"i:o:s:n:m:ecx:h")) != EOF){
		switch(opt){
			case 'i':
				inFileName = optarg;
				break;
			case 'o':
				outFileName = optarg;
				break;
			case 's':
				skipNum = atoi(optarg);
				break;
			case 'n':
				maxNum = atoi(optarg);
				break;
			case 'm':
				maxGbp = atoi(optarg);
				break;
			case 'e':
				mode += 1;
				break;
			case 'c':
				mode += 2;
				break;
			case 'x':
				X = atoi(optarg);
				break;
			case 'h':
			case '?':
			default:
				printHelp();
				return false;
		}
	}
	if(inFileName == "" || outFileName == "" || mode > 2 ){
		printHelp();
		return false;
	}
	return true;
}

void printHelp(){
	cerr << "\t***** " << progName << " *****\n\t- Andrew Spriggs, CSIRO Ag&Food, 2018 -\n\n";
	cerr << "Usage:\t" << progName << " -i seqFile -o outFile [options]\n\n";
	cerr << "Extracts a subset of a set of sequences, given rules.\n";
	cerr << "Options:\n";
	cerr << "\t-i seqFile\tFasta or fastq sequence file. May be .gz compressed.\n";
	cerr << "\t-o outFile\tOutput file, retains input format.\n";
	cerr << "\t-s N\t\tNumber- Skip over this first N sequences in input.\n";
	cerr << "\t-n maxNum\tNumber- Write this maximum count of sequences to output file.\n";
	cerr << "\t-m maxGbp\tNumber- Write this maximum Gbp of sequence to output file.\n";
	cerr << "\t-e\t\tExtract every Xth seq (-x below). For retaining <50% of input.\n";
	cerr << "\t-c\t\tExclude every Xth seq (-x below). For retaining >50% of input.\n";
	cerr << "\t-x X\t\tNumber- For the extract and eclude modes.\n";
	cerr << "Extract-every-X and exclude-every-X modes and mutually exclusive.\n\n";
}

