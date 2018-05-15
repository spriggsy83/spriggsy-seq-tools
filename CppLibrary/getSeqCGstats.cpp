#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include "SeqReader.h"
using namespace std;

/* By Andrew Spriggs, CSIRO Ag&Food, 2018 */
/* andrew.spriggs@csiro.au */
/* https://github.com/spriggsy83 */

const char progName[] = "getSeqCGstats";

bool getInputs(int argc, char* argv[], string& inFileName, string& outFileName);

int main(int argc,char *argv[]){

	string inFileName;
	string outFileName;
	
	if(!getInputs(argc, argv, inFileName, outFileName)){
		cerr << "Process aborted.\n";
		return 0;
	}
	
	ofstream outfile(outFileName.c_str());
	if(!outfile.is_open()){
		cerr << "Unable to open output file " << outFileName << "!\nProcess aborted.\n";
		return 0;
	}
	outfile << "ID\tLength\tCG%\tA\tT\tC\tG\n";
	
	SeqReader inFile(inFileName);
	
	while(inFile.nextSeq()){
		
		outfile << inFile.getSeqID() << "\t" << inFile.getSeqLen() << "\t";
		
		string seq = inFile.getSeq();
		
		int aCount = 0;
		int tCount = 0;
		int cCount = 0;
		int gCount = 0;
		
		for(int i=0; i<seq.length(); i++){
			switch (seq[i]){
				case 'A':
				case 'a':
					aCount++;
					break;
				case 'T':
				case 't':
					tCount++;
					break;
				case 'C':
				case 'c':
					cCount++;
					break;
				case 'G':
				case 'g':
					gCount++;
			}
		}
		
		double cgPC = (cCount+gCount) / (double)seq.length() * 100.00;
		outfile.setf(ios::fixed);
		outfile << setprecision(2) << cgPC << "\t";
		outfile << aCount << "\t";
		outfile << tCount << "\t";
		outfile << cCount << "\t";
		outfile << gCount << "\n";
	}

	return 0;
}

bool getInputs(int argc, char* argv[], string& inFileName, string& outFileName){
	if(argc != 3){
		cerr << "\t***** " << progName << " *****\n\t- Andrew Spriggs, CSIRO Ag&Food, 2018 -\n";
		cerr << "Profiles sequences for GC% and base counts\n";
		cerr << "Input files may be fasta or fastq and may be .gz compressed.\n";
		cerr << "Command line usage:\n" << argv[0] << " <in file> <out file>\n";
		return false;
	}
	
	inFileName = argv[1];
	outFileName = argv[2];

	return true;
}
