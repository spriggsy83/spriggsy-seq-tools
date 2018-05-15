#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include "SeqReader.h"
using namespace std;

/* By Andrew Spriggs, CSIRO Ag&Food, 2018 */
/* andrew.spriggs@csiro.au */
/* https://github.com/spriggsy83 */

/** Produces a table of sequence sizes and their counts per input **/

const char progName[] = "getSeqSizeChart";

bool getInputs(int argc, char* argv[], vector<string>& inFileNames);

int main(int argc,char *argv[]){

	vector<string> inFileNames;	
	if(!getInputs(argc, argv, inFileNames)){
		cerr << "Process aborted.\n";
		return 0;
	}
	const int numSamples = inFileNames.size();
	vector<vector<int> > seqSizes(numSamples, vector<int>(1, 0));
	
	for(int fileNum = 0; fileNum < numSamples; fileNum++){
		
		SeqReader inFile(inFileNames[fileNum]);
				
		cerr << "Processing " << inFileNames[fileNum] << "\n";
		
		while(inFile.nextSeq()){
			int len = inFile.getSeqLen();
			if(len >= seqSizes[0].size()){
				for(int i = 0; i < numSamples; i++){
					for(int j = seqSizes[i].size(); j <= len; j++){
						seqSizes[i].push_back(0);
					}
				}
			}
			seqSizes[fileNum][len] = seqSizes[fileNum][len] + 1;
		}
	}
	
	cout << "Size";
	for(int fileNum = 0; fileNum < numSamples; fileNum++){
		cout << "\t" << inFileNames[fileNum];
	}
	cout << "\n";
	for(int len = 0; len < seqSizes[0].size(); len++){
		cout << len;
		for(int fileNum = 0; fileNum < numSamples; fileNum++){
			cout << "\t" << seqSizes[fileNum][len];
		}
		cout << "\n";
	}
	cout << "\n";

	return 0;
}

bool getInputs(int argc, char* argv[], vector<string>& inFileNames){
	if(argc < 2){
		cerr << "\t***** " << progName << " *****\n\t- Andrew Spriggs, CSIRO Ag&Food, 2018 -\n";
		cerr << "Produces a table of sequence sizes and their counts per input.\n";
		cerr << "Input files may be fasta or fastq and may be .gz compressed.\n";
		cerr << "Command line usage:\n" << argv[0] << " <in seq file> [more in files]\n";
		return false;
	}
	
	for(int i = 1; i < argc; i++){
		string aFileName(argv[i]);
		inFileNames.push_back(aFileName);
	}
	return true;
}

