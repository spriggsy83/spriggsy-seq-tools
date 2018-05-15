#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include "SeqReader.h"
using namespace std;

/* By Andrew Spriggs, CSIRO Ag&Food, 2018 */
/* andrew.spriggs@csiro.au */
/* https://github.com/spriggsy83 */

const char progName[] = "getSeqQCStats";

bool getInputs(int argc, char* argv[], int& profileEvery, vector<string>& inFileNames);

int main(int argc,char *argv[]){

	vector<string> inFileNames;
	int profileEvery;
	
	if(!getInputs(argc, argv, profileEvery, inFileNames)){
		cerr << "Process aborted.\n";
		return 0;
	}
	
	const int numSamples = inFileNames.size();
	vector<unsigned long> counts(numSamples, 0);
	vector<unsigned long long> totalLens(numSamples, 0);
	vector<double> avLens(numSamples, 0);
	vector<double> avGCs(numSamples, 0);
	vector<double> avQuals(numSamples, 0);
	vector<double> avMidPQuals(numSamples, 0);
	vector<double> avEndQuals(numSamples, 0);
	
	for(int fileNum = 0; fileNum < numSamples; fileNum++){
		
		SeqReader inFile(inFileNames[fileNum]);
				
		cerr << "Processing " << inFileNames[fileNum] << "\n";
		
		long double allGCs = 0;
		long double allQuals = 0;
		unsigned int allMidPQuals = 0;
		unsigned int allEndQuals = 0;

		int profileCounter = profileEvery - 1;
		unsigned long totProfiled = 0;
		
		while(inFile.nextSeq()){
			int len = inFile.getSeqLen();
			counts[fileNum]++;
			totalLens[fileNum] += len;
			profileCounter++;

			if(profileCounter == profileEvery){
				profileCounter = 0;
				totProfiled++;

				string thisSeq = inFile.getSeq();
				int thisGCraw = 0;
				for(int i = 0; i < thisSeq.length(); i++){
					if(thisSeq[i] == 'G' || thisSeq[i] == 'g' || thisSeq[i] == 'C' || thisSeq[i] == 'c'){
						thisGCraw++;
					}
				}
				allGCs += thisGCraw / (long double)thisSeq.length() * 100.0;

				string thisQual = inFile.getSeqQual();
				if(thisQual.length() > 0){
					int thisQualTot = 0;
					for(int i = 0; i < thisQual.length(); i++){
						thisQualTot += int(thisQual[i]) - 33;
					}
					allQuals += thisQualTot / (long double)thisQual.length();
					allEndQuals += int(thisQual[thisQual.length()-1]) - 33;
					allMidPQuals += int(thisQual[int(thisQual.length()/2)]) - 33;
				}else{
					allQuals += 40;
					allMidPQuals += 40;
					allEndQuals += 40;
				}
			}

		}
		
		avLens[fileNum] = totalLens[fileNum] / (long double)counts[fileNum];
		avGCs[fileNum] = allGCs / (long double)totProfiled;
		avQuals[fileNum] = allQuals / (long double)totProfiled;	
		avMidPQuals[fileNum] = allMidPQuals / (long double)totProfiled;	
		avEndQuals[fileNum] = allEndQuals / (long double)totProfiled;	
	}
	
	cout.setf(ios::fixed);
	cout << setprecision(2);

	cout << "\tNumSeqs\tAvLen\tTotLen\tAvQual\tAvMidQual\tAvEndQual\tAvGC\n";
	for(int i = 0; i < numSamples; i++){
		cout << inFileNames[i];
		cout << "\t" << counts[i];
		cout << "\t" << avLens[i];
		cout << "\t" << totalLens[i];
		cout << "\t" << avQuals[i];
		cout << "\t" << avMidPQuals[i];
		cout << "\t" << avEndQuals[i];
		cout << "\t" << avGCs[i] << "%";
		cout << "\n";
	}
	return 0;
}

bool getInputs(int argc, char* argv[], int& profileEvery, vector<string>& inFileNames){
	if(argc < 3){
		cerr << "\t***** " << progName << " *****\n\t- Andrew Spriggs, CSIRO Ag&Food, 2018 -\n";
		cerr << "Profiles sequences for count, total/average/median bp length, av. GC%, fq phred scores...\n";
		cerr << "Input files may be fasta or fastq and may be .gz compressed.\n";
		cerr << "Command line usage:\n" << argv[0] << " <profile every X seqs> <in file> [more in files]\n";
		return false;
	}
	
	profileEvery = atoi(argv[1]);
	for(int i = 2; i < argc; i++){
		string aFileName(argv[i]);
		inFileNames.push_back(aFileName);
	}
	return true;
}

