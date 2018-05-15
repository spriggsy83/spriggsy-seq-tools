#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <set>
#include <vector>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include "SeqReader.h"
using namespace std;

/* By Andrew Spriggs, CSIRO Ag&Food, 2018 */
/* andrew.spriggs@csiro.au */
/* https://github.com/spriggsy83 */

const char progName[] = "excludeSeqsBySAM";
/** Filters a fasta/fastq file to exclude sequences not listed as aligned in a SAM file. **/

bool getInputs(int argc, char* argv[], string& inSeqsFileName, string& inSAMFileName, string& outFileName);

int main(int argc,char *argv[]){

	string inSeqsFileName;
	string inSAMFileName;
	string outFileName;
	set<string> readIDs;

	if(!getInputs(argc, argv, inSeqsFileName, inSAMFileName, outFileName)){
		cerr << "Process aborted.\n";
		return 1;
	}
	
	// Prepare output 
	ofstream outfile(outFileName.c_str());
	if(!outfile.is_open()){
		cerr << "Unable to open output file " << outFileName << "!\nProcess aborted.\n";
		return 1;
	}
	
	// Read SAM file into readIDs list
	ifstream fileifs;
	boost::iostreams::filtering_istream infile;
	bool gzipFile = false;
	if(inSAMFileName.find("gz", inSAMFileName.length()-3) != string::npos || 
			inSAMFileName.find("GZ", inSAMFileName.length()-3) != string::npos){
		gzipFile = true;
	}
	if(gzipFile){
		fileifs.open(inSAMFileName.c_str(), ios_base::in | ios_base::binary);
	}else{
		fileifs.open(inSAMFileName.c_str(), ios_base::in );
	}
	if (!fileifs.is_open()){
		cerr << "Unable to open SAM file " << inSAMFileName << "!\n";
		outfile.close();
		return 1;
	}else if (!fileifs.good()){
		cerr << "SAM file " << inSAMFileName << " is empty!\n";
		outfile.close();
		fileifs.close();
		return 1;
	}
	
	try {
		if(gzipFile){
			infile.push(boost::iostreams::gzip_decompressor());
		}
		infile.push(fileifs);

		string line;
		while(getline(infile, line)){
			stringstream linestream(line);
			vector<string> lineParts;
			lineParts.reserve(11);
			string aLinePart;
			// Tab separated split
			while(getline(linestream, aLinePart, '\t')){
				lineParts.push_back(aLinePart);
			}
			if(lineParts.size() == 11){
				// readID == [0], refID == [2], start == [3], cigar == [5], readSeq == [9]
				readIDs.insert(lineParts[0]);
			}
		}
	}
	catch(const boost::iostreams::gzip_error& e) {
		cerr << "Error while reading SAM file " << inSAMFileName << endl;
		cerr << e.what() << endl;
	}
	fileifs.close();
	infile.reset();

	cout << "Read " << readIDs.size() << " unique read IDs to exclude.\n";
	
	// Read sequence file and write filtered result

	SeqReader inFile(inSeqsFileName);
	unsigned long seqCount = 0;
	unsigned long printed = 0;
	while(inFile.nextSeq()){
		seqCount++;
		string ID = inFile.getSeqID();
		int notSpacePos = ID.find_first_not_of(" \t");
		if(notSpacePos != string::npos){
			if(notSpacePos > 0){
				ID = ID.substr(notSpacePos);
			}
			int spacePos = ID.find_first_of(" \t");
			if(spacePos != string::npos){
				ID = ID.substr(0, spacePos);
			}
			if(readIDs.count(ID) == 0){
				printed++;
				outfile << inFile.toString();
			}
		}
	}
	
	outfile.close();
	
	cout << "Read " << seqCount << " reads.\n";
	cout << "Retained " << printed << " reads.\n";
	
	return 0;
}

bool getInputs(int argc, char* argv[], string& inSeqsFileName, string& inSAMFileName, string& outFileName){
	if(argc != 4){
		cerr << "\t***** " << progName << " *****\n\t- Andrew Spriggs, CSIRO Ag&Food, 2018 -\n";
		cerr << "Filters a fasta/fastq file to exclude sequences not listed as aligned in a SAM file.\n";
		cerr << "Input files may be .gz compressed.\n";
		cerr << "Warning: Sequence IDs will be truncated at a space character.\n";
		cerr << "Correct command line usage:\n" << argv[0] << " <in seqs file> <in SAM file> <out seqs file>\n";
		return false;
	}
	
	inSeqsFileName = argv[1];
	inSAMFileName = argv[2];
	outFileName = argv[3];
	return true;
}

