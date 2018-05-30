#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <utility>
#include <cstdlib>
#include <boost/regex.hpp>
#include "SeqReader.h"
using namespace std;

/* By Andrew Spriggs, CSIRO Ag&Food, 2018 */
/* andrew.spriggs@csiro.au */
/* https://github.com/spriggsy83 */

const char progName[] = "getSubSeqs";

bool getInputs(int argc, char* argv[], string& inFileName, string& inCoordsFileName, string& outFileName);

int main(int argc,char *argv[]){

	string inFileName;
	string inCoordsFileName;
	string outFileName;
	map< string, vector< pair<int,int> > > coordList;
	map< string, bool > printedSeqs;
	
	if(!getInputs(argc, argv, inFileName, inCoordsFileName, outFileName)){
		cerr << "Process aborted.\n";
		return 0;
	}


		/** Coordinate collection **/
	int coordCount = 0;
	ifstream coordfile;
	string line;
	coordfile.open(inCoordsFileName.c_str());
	if (!coordfile.is_open()){
		cerr << "Unable to open coordinates file " << inCoordsFileName << "!\n";
		return false;
	}else if (!coordfile.good()){
		cerr << "Coordinates file " << inCoordsFileName << " is empty!\n";
		coordfile.close();
		return false;
	}
	string prevSeqID ("");
	while(getline(coordfile, line)){
		boost::regex pattern1 ("^(\\S+):(\\d+)-(\\d+).*");
		boost::regex pattern2 ("^(\\S+)\t(\\d+)\t(\\d+).*");
		boost::regex pattern3 ("^(\\d+)-(\\d+)$");
		boost::regex pattern4 ("^>(\\S+).*");
		boost::smatch matchedPats;
		if(boost::regex_match(line, matchedPats, pattern1)){
			string seqID = matchedPats[1].str();
			pair<int,int> coords ( atoi(matchedPats[2].str().c_str()), atoi(matchedPats[3].str().c_str()) );
			if(coordList.count(seqID) == 0){
				coordList[seqID] = vector< pair<int,int> >();
			}
			coordList[seqID].push_back( coords );
			coordCount++;
			printedSeqs[seqID] = false;

		}else if(boost::regex_match(line, matchedPats, pattern2)){
			string seqID = matchedPats[1].str();
			pair<int,int> coords ( atoi(matchedPats[2].str().c_str()), atoi(matchedPats[3].str().c_str()) );
			if(coordList.count(seqID) == 0){
				coordList[seqID] = vector< pair<int,int> >();
			}
			coordList[seqID].push_back( coords );
			coordCount++;
			printedSeqs[seqID] = false;

		}else if(boost::regex_match(line, matchedPats, pattern3)){
			if(prevSeqID.length() > 0){
				pair<int,int> coords ( atoi(matchedPats[1].str().c_str()), atoi(matchedPats[2].str().c_str()) );
				if(coordList.count(prevSeqID) == 0){
					coordList[prevSeqID] = vector< pair<int,int> >();
				}
				coordList[prevSeqID].push_back( coords );
				coordCount++;
				printedSeqs[prevSeqID] = false;
			}

		}else if(boost::regex_match(line, matchedPats, pattern4)){
			prevSeqID = matchedPats[1].str();
		}
	}
	coordfile.close();
	cout << "Loaded " << coordCount << " coordinates on " << coordList.size() << " sequences.\n";
	

		/** Result writing **/
	ofstream outfile(outFileName.c_str());
	if(!outfile.is_open()){
		cerr << "Unable to open output file " << outFileName << "!\nProcess aborted.\n";
		return 0;
	}
	
	SeqReader inSeqs(inFileName);
	
	int readcount = 0;
	int writecount = 0;
	while(inSeqs.nextSeq()){
		readcount++;
		string seqID = inSeqs.getSeqID();
		int seqLen = inSeqs.getSeqLen();

		if(coordList.count(seqID) > 0){
			printedSeqs[seqID] = true;

			for(vector< pair<int,int> >::iterator aCoord=coordList[seqID].begin(); aCoord!=coordList[seqID].end(); ++aCoord){
				int start = aCoord->first;
				int end = aCoord->second;
				if(start < 1 || start > seqLen || end < 1 || end > seqLen ){
					cout << "Invalid coords, " << start << "-" << end << ", on " << seqID << ", length " << seqLen << "\n"; 
				}else{
					outfile << inSeqs.toStringSubseq(start, end);
					writecount++;
				}
			}
		}
	}
	outfile.close();

	for(map<string,bool>::iterator wasItPrinted=printedSeqs.begin(); wasItPrinted!=printedSeqs.end(); ++wasItPrinted){
		if(wasItPrinted->second == false){
			cout << "ID not found in sequence file: " << wasItPrinted->first << "\n";
		}
	}
	
	cout << "Read " << readcount << " sequences.\n";
	cout << "Wrote " << writecount << " sub-sequences.\n";
	return 0;
}

bool getInputs(int argc, char* argv[], string& inFileName, string& inCoordsFileName, string& outFileName){
	if(argc != 4){
		cerr << "\t***** " << progName << " *****\n\t- Andrew Spriggs, CSIRO Ag&Food, 2018 -\n";
		cerr << "Will output sub-sequences, from a list of sequence coordinate ranges.\n";
		cerr << "Input sequence file may be fasta or fastq and may be .gz compressed.\n\n";
		cerr << "Coordinates for sub-sequence ranges may be in various forms:\n";
		cerr << "  SeqID:start-end\n";
		cerr << "  SeqID\tstart\tend\t(tab-separated)\n";
		cerr << "  >SeqID\n  start-end\n";
		cerr << "Coordinate ranges are inclusive and count from _1_\n\n";
		cerr << "Command line usage:\n" << argv[0] << " <in seq file> <in coords file> <outfile>\n";
		return false;
	}
	
	inFileName = argv[1];
	inCoordsFileName = argv[2];
	outFileName = argv[3];
	return true;
}

