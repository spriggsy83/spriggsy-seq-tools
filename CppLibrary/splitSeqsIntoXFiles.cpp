#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <vector>
#include <sstream>
#include "SeqReader.h"
using namespace std;

/* By Andrew Spriggs, CSIRO Ag&Food, 2018 */
/* andrew.spriggs@csiro.au */
/* https://github.com/spriggsy83 */

const char progName[] = "splitSeqsIntoXFiles";

bool getInputs(int argc, char* argv[], string& inFileName, string& outFilePref, int& numFiles);

int main(int argc,char *argv[]){

	string inFileName;
	string outFilePref;
	int numFiles = 1;
	vector<ofstream*> outFiles;
	
	if(!getInputs(argc, argv, inFileName, outFilePref, numFiles)){
		cerr << "Process aborted.\n";
		return 1;
	}
	SeqReader inFile(inFileName);
	
	string outFileExt = inFile.fileModeString();
	outFiles.resize(numFiles);
	for(int i=0; i < numFiles; i++){
		stringstream outFileNameS;
		outFileNameS << outFilePref << "-pt" << (i+1) << "." << outFileExt;
		outFiles[i] = new ofstream(outFileNameS.str().c_str());
		if(!outFiles[i]->is_open()){
			cerr << "Unable to open output file " << outFileNameS.str() << "!\nProcess aborted.\n";
			for(int j=0; j < i; j++){
				outFiles[j]->close();
				delete outFiles[j];
			}
			return 1;
		}
		cout << "Opened output: " << outFileNameS.str() << "\n";
	}
	
	int counter = 0;
	unsigned long numPrinted = 0;
	while(inFile.nextSeq()){
		*outFiles[counter] << inFile.toString();
		(*outFiles[counter]).flush();
		numPrinted++;
		counter++;
		if(counter == numFiles){
			counter = 0;
		}
	}
	
	for(int i=0; i < outFiles.size(); i++){
		outFiles[i]->close();
		delete outFiles[i];
	}
	
	cout << "Num printed\t" << numPrinted << "\n";
	
	return 0;
}

bool getInputs(int argc, char* argv[], string& inFileName, string& outFileName, int& selectNum){
	if(argc != 4){
		cerr << "\t***** " << progName << " *****\n\t- Andrew Spriggs, CSIRO Ag&Food, 2018 -\n";
		cerr << "Command line usage:\n" << argv[0] << " <infile> <outfile-prefix> <numFiles>\n";
		cerr << "Will divide a sequence file into [numFiles] pieces.\n";
		cerr << "Out file names = [outfile-prefix]-ptX.[infile's-suffix]\n";
		cerr << "Input files may be fasta or fastq and may be .gz compressed.\n";
		return false;
	}
	
	inFileName = argv[1];
	outFileName = argv[2];
	selectNum = atoi(argv[3]);
	return true;
}

