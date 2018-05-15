#include <iostream>
#include <fstream>
#include <string>
#include "SeqReader.h"
using namespace std;

/* By Andrew Spriggs, CSIRO Ag&Food, 2018 */
/* andrew.spriggs@csiro.au */
/* https://github.com/spriggsy83 */

const char progName[] = "getSeqSizeList";

bool getInputs(int argc, char* argv[], string& inFileName);

int main(int argc,char *argv[]){

	string inFileName;
	
	if(!getInputs(argc, argv, inFileName)){
		cerr << "Process aborted.\n";
		return 0;
	}
	SeqReader inFile(inFileName);
				
	while(inFile.nextSeq()){
		int len = inFile.getSeqLen();
		cout << inFile.getSeqID() << "\t" << inFile.getSeqLen() << "\n";
	}

	return 0;
}

bool getInputs(int argc, char* argv[], string& inFileName){
	if(argc < 1){
		cerr << "\t***** " << progName << " *****\n\t- Andrew Spriggs, CSIRO Ag&Food, 2018 -\n";
		cerr << "Print list of sequence bp lengths to stdout\n";
		cerr << "Input files may be fasta or fastq and may be .gz compressed.\n";
		cerr << "Command line usage:\n" << argv[0] << " <in seq file>\n";
		return false;
	}
	string aFileName(argv[1]);
	inFileName = aFileName;
	return true;
}

