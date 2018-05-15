#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include "SeqReader.h"
using namespace std;

/* By Andrew Spriggs, CSIRO Ag&Food, 2018 */
/* andrew.spriggs@csiro.au */
/* https://github.com/spriggsy83 */

const char progName[] = "extractEveryXSeq";

bool getInputs(int argc, char* argv[], string& inFileName, string& outFileName, int& selectNum);

int main(int argc,char *argv[]){

	string inFileName;
	string outFileName;
	int selectNum = 1;
	
	unsigned long numPrinted = 0;
	unsigned long numRejected = 0;
	
	
	if(!getInputs(argc, argv, inFileName, outFileName, selectNum)){
		cerr << "Process aborted.\n";
		return 0;
	}
	
	ofstream outfile(outFileName.c_str());
	if(!outfile.is_open()){
		cerr << "Unable to open output file " << outFileName << "!\nProcess aborted.\n";
		return 0;
	}
	
	SeqReader inFile(inFileName);
	
	int counter = 0;
	while(inFile.nextSeq()){
		counter++;
		if(counter == selectNum){
			outfile << inFile.toString();
			numPrinted++;
			counter = 0;
		}else{
			numRejected++;
		}
	}
	
	cout << "Num printed\t" << numPrinted << "\n";
	cout << "Num rejected\t" << numRejected << "\n";
	return 0;
}

bool getInputs(int argc, char* argv[], string& inFileName, string& outFileName, int& selectNum){
	if(argc != 4){
		cerr << "\t***** " << progName << " *****\n\t- Andrew Spriggs, CSIRO Ag&Food, 2018 -\n";
		cerr << "Will output a subset of a sequence file, retaining every Xth sequence.\n";
		cerr << "Input files may be fasta or fastq and may be .gz compressed.\n";
		cerr << "Command line usage:\n" << argv[0] << " <in seq file <outfile> <X>\n";
		return false;
	}
	
	inFileName = argv[1];
	outFileName = argv[2];
	selectNum = atoi(argv[3]);
	return true;
}

