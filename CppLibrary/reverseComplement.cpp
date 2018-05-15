#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include "SeqReader.h"
using namespace std;

/* By Andrew Spriggs, CSIRO Ag&Food, 2018 */
/* andrew.spriggs@csiro.au */
/* https://github.com/spriggsy83 */

/** Reverse complements each sequence from a file > FASTA */

const char progName[] = "reverseComplement";

bool getInputs(int argc, char* argv[], string& inFileName, string& outFileName);

int main(int argc,char *argv[]){

	string inFileName;
	string outFileName;
	
	if(!getInputs(argc, argv, inFileName, outFileName)){
		cerr << "Process aborted.\n";
		return 0;
	}
	
	ofstream outFile(outFileName.c_str());
	if(!outFile.is_open()){
		cerr << "Unable to open output file " << outFileName << "!\nProcess aborted.\n";
		return 0;
	}
	
	SeqReader inFile(inFileName);
	
	while(inFile.nextSeq()){
		outFile << ">" << inFile.getSeqID() << "\n";
		outFile << inFile.revComp() << "\n";
	}
	
	outFile.close();
	return 0;
}

bool getInputs(int argc, char* argv[], string& inFileName, string& outFileName){
	if(argc != 3){
		cerr << "\t***** " << progName << " *****\n\t- Andrew Spriggs, CSIRO Ag&Food, 2018 -\n";
		cerr << "Will output a reverse complemented set of sequences in FASTA format.\n";
		cerr << "Input file may be fasta or fastq and may be .gz compressed.\n";
		cerr << "Command line usage:\n" << argv[0] << " <infile> <out fasta file>\n";
		return false;
	}
	
	inFileName = argv[1];
	outFileName = argv[2];
	return true;
}
