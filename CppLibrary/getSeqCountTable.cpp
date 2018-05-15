#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include "SeqReader.h"
using namespace std;

/* By Andrew Spriggs, CSIRO Ag&Food, 2018 */
/* andrew.spriggs@csiro.au */
/* https://github.com/spriggsy83 */

const char progName[] = "getSeqCountTable";

bool getInputs(int argc, char* argv[], vector<string>& inFileNames, string& outFileName, bool& filterPoly, bool& singletons);
bool isPolySeq(const string& seq);

int main(int argc,char *argv[]){

	vector<string> inFileNames;
	string outFileName;
	bool filterPoly;
	bool singletons;
	
	if(!getInputs(argc, argv, inFileNames, outFileName, filterPoly, singletons)){
		cerr << "Process aborted.\n";
		return 0;
	}
	
	const int numSamples = inFileNames.size();
	map<string, vector<unsigned long> > counts;
	
	for(int fileNum = 0; fileNum < numSamples; fileNum++){
		
		SeqReader inFile(inFileNames[fileNum]);
				
		cout << "Processing " << inFileNames[fileNum] << "\n";
		
		while(inFile.nextSeq()){
			string seq = inFile.getSeq();
			
			map<string, vector<unsigned long> >::iterator seqRec = counts.find(seq);
			
			if(seqRec != counts.end()){
				(*seqRec).second[fileNum] = (*seqRec).second[fileNum] + 1;
				(*seqRec).second[numSamples] = (*seqRec).second[numSamples] + 1;
			}else{
				vector<unsigned long> countVec(numSamples+1, 0);
				countVec[fileNum] = 1;
				countVec[numSamples] = 1;
				counts[seq] = countVec;
			}
		}
	}
	
	cout << "Writing output...\n";
	if(filterPoly){
		cout << "...and filtering out poly-N sequences...\n";
	}
	
	ofstream outfile(outFileName.c_str());
	if(!outfile.is_open()){
		cerr << "Unable to open output file " << outFileName << "!\nProcess aborted.\n";
		return 0;
	}
	
	unsigned long seqNum = 1;
	outfile << "ID\tSeq";
	for(int fileNum = 0; fileNum < numSamples; fileNum++){
		outfile << "\t" << inFileNames[fileNum];
	}
	outfile << "\n";
	
	for (map<string, vector<unsigned long> >::iterator seqRec = counts.begin(); seqRec!=counts.end(); ++seqRec){
		
		if( singletons || (*seqRec).second[numSamples] > 1 ){
			if( !filterPoly || !isPolySeq((*seqRec).first) ){
				outfile << seqNum << "\t" << (*seqRec).first;
				for(int fileNum = 0; fileNum < numSamples; fileNum++){
					outfile << "\t" << (*seqRec).second[fileNum];
				}
				outfile << "\n";
			}
		}
		seqNum++;
	}
	
	outfile.close();

	return 0;
}

bool isPolySeq(const string& seq){
	
	const float maxSingle = seq.length() * 2.0 / 3.0;
	const float maxDuos = seq.length() * 2.0 / 6.0;
	set<string> tested;
	
	for(int i=0; i<seq.length(); i++){
		string base = (seq.substr(i, 1));
		if(tested.find(base) == tested.end()){
			tested.insert(base);
			int count = 1;
			for(int j=i+1; j<seq.length(); j++){
				if(seq[j] == seq[i]){
					count++;
				}
			}
			if(count >= maxSingle){
				return true;
			}
		}
		if(i>0){
			string duo = (seq.substr(i-1, 2));
			if(tested.find(duo) == tested.end()){
				tested.insert(duo);
				int count = 1;
				for(int j=i+2; j<seq.length(); j=j+2){
					if(seq[j-1] == seq[i-1] && seq[j] == seq[i]){
						count++;
					}
				}
				if(count >= maxDuos){
					return true;
				}
			}
		}
	}
	
	return false;
}

bool getInputs(int argc, char* argv[], vector<string>& inFileNames, string& outFileName, bool& filterPoly, bool& singletons){
	if(argc < 5){
		cerr << "\t***** " << progName << " *****\n\t- Andrew Spriggs, CSIRO Ag&Food, 2018 -\n";
		cerr << "Produces a table of occurances of individual sequences per input.\n";
		cerr << "Input files may be fasta or fastq and may be .gz compressed.\n";
		cerr << "Command line usage:\n" << argv[0] << " <out tab file> <filter poly-N T/F> <remove singletons T/F> <in seq file> [more in files]\n";
		return false;
	}
	
	outFileName = argv[1];
	string filterAns = argv[2];
	if(filterAns[0] == 'T' || filterAns[0] == 't'){
		filterPoly = true;
	}else if(filterAns[0] == 'F' || filterAns[0] == 'f'){
		filterPoly = false;
	}else{
		cerr << "\t***** " << progName << " *****\n\t- Andrew Spriggs, CSIRO Ag&Food, 2018 -\n";
		cerr << "Correct command line usage:\n" << argv[0] << " <out tab file> <filter poly-N T/F> <remove singletons T/F> <in seq file> [more in files]\n";
		return false;
	}
	string singlesAns = argv[3];
	if(singlesAns[0] == 'T' || singlesAns[0] == 't'){
		singletons = false;
	}else if(singlesAns[0] == 'F' || singlesAns[0] == 'f'){
		singletons = true;
	}else{
		cerr << "\t***** " << progName << " *****\n\t- Andrew Spriggs, CSIRO Ag&Food, 2018 -\n";
		cerr << "Correct command line usage:\n" << argv[0] << " <out tab file> <filter poly-N T/F> <remove singletons T/F> <in seq file> [more in files]\n";
		return false;
	}
	
	for(int i = 4; i < argc; i++){
		string aFileName(argv[i]);
		inFileNames.push_back(aFileName);
	}
	return true;
}

