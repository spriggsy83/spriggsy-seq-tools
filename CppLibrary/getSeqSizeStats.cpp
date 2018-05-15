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

const char progName[] = "getSeqSizeStats";

bool getInputs(int argc, char* argv[], vector<string>& inFileNames);

int main(int argc,char *argv[]){

	vector<string> inFileNames;
	
	if(!getInputs(argc, argv, inFileNames)){
		cerr << "Process aborted.\n";
		return 0;
	}
	
	const int numSamples = inFileNames.size();
	vector<unsigned long> counts(numSamples, 0);
	vector<unsigned long long> totals(numSamples, 0);
	vector<unsigned int> mins(numSamples, -1);
	vector<unsigned int> maxs(numSamples, -1);
	string minId = "";
	string maxId = "";
	vector<double> avs(numSamples, 0);
	vector<unsigned int> medians(numSamples, 0);
	vector<unsigned int> n50s(numSamples, 0);
	
	for(int fileNum = 0; fileNum < numSamples; fileNum++){
		
		SeqReader inFile(inFileNames[fileNum]);
				
		cerr << "Processing " << inFileNames[fileNum] << "\n";
		
		vector<int> allLengths;
		
		while(inFile.nextSeq()){
			int len = inFile.getSeqLen();
			allLengths.push_back(len);
			counts[fileNum]++;
			totals[fileNum] += len;
			if(len < mins[fileNum] || mins[fileNum] == -1){
				mins[fileNum] = len;
				minId = inFile.getSeqID();
			}
			if(len > maxs[fileNum] || maxs[fileNum] == -1){
				maxs[fileNum] = len;
				maxId = inFile.getSeqID();
			}
		}
		
		avs[fileNum] = totals[fileNum] / (long double)counts[fileNum];
		
		sort (allLengths.begin(), allLengths.end());
		if(counts[fileNum] % 2 == 1){
			medians[fileNum] = allLengths[int(counts[fileNum]/2)];
		}else{
			medians[fileNum] = (allLengths[counts[fileNum]/2] + allLengths[counts[fileNum]/2 - 1]) / 2;
		}

		unsigned long long n50Tot = 0;
		for(unsigned long i = 0; i < allLengths.size() && n50Tot < (totals[fileNum] / 2.0); i++){
			n50Tot += allLengths[i];
			n50s[fileNum] = allLengths[i];
		}
	}
	
	for(int i = 0; i < numSamples; i++){
		cout << "\t" << inFileNames[i];
	}
	cout << "\nNumber of sequences";
	for(int i = 0; i < numSamples; i++){
		cout << "\t" << counts[i];
	}
	cout << "\nCombined length";
	for(int i = 0; i < numSamples; i++){
		cout << "\t" << totals[i];
	}
	cout << "\nMinimum sequence length";
	for(int i = 0; i < numSamples; i++){
		cout << "\t" << mins[i];
	}
	if(numSamples == 1){
		cout << "\t(" << minId << ")";
	}
	cout << "\nAverage sequence length";
	for(int i = 0; i < numSamples; i++){
		cout << "\t" << avs[i];
	}
	cout << "\nMedian sequence length";
	for(int i = 0; i < numSamples; i++){
		cout << "\t" << medians[i];
	}
	cout << "\nN50 sequence length";
	for(int i = 0; i < numSamples; i++){
		cout << "\t" << n50s[i];
	}
	cout << "\nMaximum sequence length";
	for(int i = 0; i < numSamples; i++){
		cout << "\t" << maxs[i];
	}
	if(numSamples == 1){
		cout << "\t(" << maxId << ")";
	}
	cout << "\n";

	return 0;
}

bool getInputs(int argc, char* argv[], vector<string>& inFileNames){
	if(argc < 2){
		cerr << "\t***** " << progName << " *****\n\t- Andrew Spriggs, CSIRO Ag&Food, 2018 -\n";
		cerr << "Profiles sequences for total/average/median/min/max bp lengths.\n";
		cerr << "Input files may be fasta or fastq and may be .gz compressed.\n";
		cerr << "Command line usage:\n" << argv[0] << " <in file> [more in files]\n";
		return false;
	}
	
	for(int i = 1; i < argc; i++){
		string aFileName(argv[i]);
		inFileNames.push_back(aFileName);
	}
	return true;
}

