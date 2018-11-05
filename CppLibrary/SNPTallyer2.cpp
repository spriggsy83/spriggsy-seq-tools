#include <omp.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <string>
#include <cstring>
#include <ctype.h>
#include <sstream>
#include <algorithm>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include "SNPTallyer2.h"

/* By Andrew Spriggs, CSIRO Ag&Food, 2018 */
/* andrew.spriggs@csiro.au */
/* https://github.com/spriggsy83 */

/*** Initialise with no defaults 
**/
SNPTallyer::SNPTallyer(const vector<string>& aLabelsList, 
						const vector<string>& aSAMFileNamesList, 
						const vector<string>& aSNPFileNamesList, 
						const string& aRefSeqFileName, 
						const string& aOutTabFileName, 
						const int aOutFormat, 
						const int aReadDepthMin, 
						const int aEdgeBuffer){
	prepareSNPTallyer(aLabelsList, 
					aSAMFileNamesList, 
					aSNPFileNamesList, 
					aRefSeqFileName, 
					aOutTabFileName, 
					aOutFormat, 
					aReadDepthMin, 
					aEdgeBuffer);
	return;
}

/*** Actual constructor 
**/
void SNPTallyer::prepareSNPTallyer(const vector<string>& aLabelsList, 
								const vector<string>& aSAMFileNamesList, 
								const vector<string>& aSNPFileNamesList, 
								const string& aRefSeqFileName, 
								const string& aOutTabFileName, 
								const int aOutFormat, 
								const int aReadDepthMin, 
								const int aEdgeBuffer){

	filesReady = false;
	snpsPreLoaded = false;
	labels = aLabelsList;
	inSAMFileNames = aSAMFileNamesList;
	inSNPFileNames = aSNPFileNamesList;
	numSamples = labels.size();
	inRefSeqFileName = aRefSeqFileName;
	readsLoadedFor = "noneyet";
	readDepthMin = aReadDepthMin;
	edgeBuffer = aEdgeBuffer;
	outFormat = aOutFormat;
	if(outFormat < 1 || outFormat > 3){
		cerr << "Invalid output format option!\nNo SNP detection will follow.\n";
		return;
	}
	
	for(int sNum=0; sNum < numSamples; sNum++){
		omp_lock_t writelock;
		ompWriteLocks.push_back(writelock);
		omp_init_lock(&ompWriteLocks[sNum]);
	}
	
	outtabfile.open(aOutTabFileName.c_str());
	if(!outtabfile.is_open()){
		cerr << "Unable to open output file " << aOutTabFileName << "!\nNo SNP detection will follow.\n";
		return;
	}
	
	switch(outFormat){
		case 1:
			outtabfile << "RefID\tSNPCoord\tRefBase\tSNPBase";
			for(int sNum=0; sNum < numSamples; sNum++){
				outtabfile << "\t" << labels[sNum] << ".snpRds\t" << labels[sNum] << ".otherRds";
			}
			break;
		case 2:
			outtabfile << "RefID\tSNPCoord\tRowAllele";
			for(int sNum=0; sNum < numSamples; sNum++){
				outtabfile << "\t" << labels[sNum];
			}
			break;
		case 3:
			outtabfile << "RefID\tSNPCoord\tRefBase";
			for(int sNum=0; sNum < numSamples; sNum++){
				outtabfile << "\t" << labels[sNum] << ".A";
				outtabfile << "\t" << labels[sNum] << ".T";
				outtabfile << "\t" << labels[sNum] << ".C";
				outtabfile << "\t" << labels[sNum] << ".G";
			}
			break;
	}
	outtabfile << "\n";
	
	filesReady = true;
	return;
}


SNPTallyer::~SNPTallyer(){
	if(outtabfile.is_open()){
		outtabfile.close();
	}
	for(int sNum=0; sNum < numSamples; sNum++){
		omp_destroy_lock(&ompWriteLocks[sNum]);
	}
}

/*** Launch SNP tally across all reference sequences
**/
bool SNPTallyer::tallySNPs(){
	
	if(!filesReady){
		return false;
	}
	
	if(!loadSNPLists()){
		cerr << "Failed to load any starting SNPs from biokanga-align SNP files." << endl;
		return false;
	}else{
		// For each reference sequence
		SeqReader refSeqReader(inRefSeqFileName);
		while(refSeqReader.nextSeq()){
			if(! tallySNPsOnRef(refSeqReader.getSeq(), refSeqReader.getSeqID()) ){
				return false;
			}
		}
	}
	return true;
}

/*** Read Biokanga-Align SNP lists to form starting list of SNP locations
**/
bool SNPTallyer::loadSNPLists(){
	for(int sNum=0; sNum < numSamples; sNum++){

		ifstream fileifs;
		boost::iostreams::filtering_istream infile;
		bool gzipFile = false;
		if(inSNPFileNames[sNum].find("gz", inSNPFileNames[sNum].length()-3) != string::npos || 
				inSNPFileNames[sNum].find("GZ", inSNPFileNames[sNum].length()-3) != string::npos){
			gzipFile = true;
		}
		if(gzipFile){
			fileifs.open(inSNPFileNames[sNum].c_str(), ios_base::in | ios_base::binary);
		}else{
			fileifs.open(inSNPFileNames[sNum].c_str(), ios_base::in );
		}
		if (!fileifs.is_open()){
			cerr << "Unable to open SNPs file " << inSNPFileNames[sNum] << "!\n";
			return false;
		}else if (!fileifs.good()){
			cerr << "SNPs file " << inSNPFileNames[sNum] << " is empty!\n";
			fileifs.close();
			return false;
		}
		
		try {
			if(gzipFile){
				infile.push(boost::iostreams::gzip_decompressor());
			}
			infile.push(fileifs);

			cout << "Parsing SNP list from " << inSNPFileNames[sNum] << endl;

			string line;
			while(getline(infile, line)){
				stringstream linestream(line);
				vector<string> lineParts;
				lineParts.reserve(2);
				string aLinePart;
				// Tab separated split
				while(getline(linestream, aLinePart, ',')){
					lineParts.push_back(aLinePart);
				}
				
				//1,"SNP","cottonAD","A-chr11",76193,76193,1,"+",43,0.003621,3,1,"C",0,0,0,1,0,0.035157,144,3,0,0
				if(lineParts.size() == 23){
					if(lineParts[0] != "\"SNP_ID\""){
						string refID = lineParts[3].substr(1, lineParts[3].size()-2);
						stringstream coordSS(lineParts[4]);
						unsigned int coord;
						coordSS >> coord;
						stringstream mmReadsSS(lineParts[11]);
						unsigned int mmReads;
						mmReadsSS >> mmReads;
						if(mmReads >= readDepthMin){
							if(snpPreList.count(refID) == 0){
								snpPreList[refID] = set<unsigned int>();
							}
							snpPreList[refID].insert(coord);
						}
					}
				}
			}
		}
		catch(const boost::iostreams::gzip_error& e) {
			cerr << "Error while reading SNP file " << inSNPFileNames[sNum] << endl;
			cerr << e.what() << endl;
		}
		fileifs.close();
		infile.reset();
	}
	int numRefIDs = snpPreList.size();
	unsigned int numSNPs = 0;
	for (map<string, set<unsigned int> >::iterator aRef=snpPreList.begin(); aRef!=snpPreList.end(); ++aRef){
		numSNPs += aRef->second.size();
	}
	cout << "Loaded " << numSNPs << " starting SNPs over " << numRefIDs << " reference sequences." << endl;
	
	snpsPreLoaded = true;
	
	if(numSNPs == 0){
		return false;
	}
	return true;
}
	
/*** Launch SNP tally against a single reference sequence
**/
bool SNPTallyer::tallySNPsOnRef(const string& refSeq, const string& refID){
	if(!filesReady){
		return false;
	}
	if(!snpsPreLoaded){
		if(!loadSNPLists()){
			cerr << "Failed to load any starting SNPs from biokanga-align SNP files." << endl;
			return false;
		}
	}
	
	if(snpPreList.count(refID) > 0){
	
		int refSeqLen = refSeq.length();
		cout << "Working on " << refID << " (length = " << refSeqLen << ")... \n";
		
		// Load aligned reads from SAM files
		readReadsAll(refID, refSeqLen);
		
		// For each SNP on this RefSeq
		unsigned int snpPrintCount = 0;
		for(set<unsigned int>::iterator snpCoord=snpPreList[refID].begin(); snpCoord!=snpPreList[refID].end(); ++snpCoord){
			if(testSNP(*snpCoord, refSeq[*snpCoord], refID)){
				snpPrintCount++;
			}
		}
		cout << "Output SNPs at " << snpPrintCount << " coords on " << refID << endl;
	}
		
	return true;
}

/*** Load read alignments against a reference seq, from SAM files, for all samples 
**/
void SNPTallyer::readReadsAll(const string& refID, int refSeqLen){
	reads.clear();
	for(int sNum=0; sNum < numSamples; sNum++){
		reads.push_back(vector<AlignedRead>());
	}
	
	#pragma omp parallel for
	for(int sNum=0; sNum < numSamples; sNum++){	
		int totalReads = readReadsSample(inSAMFileNames[sNum], refID, reads[sNum]);
		cout << "Loaded " << totalReads << " reads aligned to " << refID << " from " << labels[sNum] << endl;
	}
	readsLoadedFor = refID;
	return;
}

/*** Load read alignments against a reference seq, from SAM files, for a single sample, adding to vector of aligned reads
 ***/
int SNPTallyer::readReadsSample(const string& inSamFileName, string refID, vector<AlignedRead>& reads){
	int count = 0;
	ifstream fileifs;
	boost::iostreams::filtering_istream infile;
	bool gzipFile = false;
	if(inSamFileName.find("gz", inSamFileName.length()-3) != string::npos || 
			inSamFileName.find("GZ", inSamFileName.length()-3) != string::npos){
		gzipFile = true;
	}
	if(gzipFile){
		fileifs.open(inSamFileName.c_str(), ios_base::in | ios_base::binary);
	}else{
		fileifs.open(inSamFileName.c_str(), ios_base::in );
	}
	if (!fileifs.is_open()){
		cerr << "Unable to open SAM file " << inSamFileName << "!\n";
		return false;
	}else if (!fileifs.good()){
		cerr << "SAM file " << inSamFileName << " is empty!\n";
		fileifs.close();
		return false;
	}
	
	try {
		if(gzipFile){
			infile.push(boost::iostreams::gzip_decompressor());
		}
		infile.push(fileifs);

		refID = "\t" + refID + "\t";
		
		string line;
		while(getline(infile, line)){
			// If refID in line
			if(line.find(refID) != string::npos){
				stringstream linestream(line);
				vector<string> lineParts;
				lineParts.reserve(11);
				string aLinePart;
				// Tab separated split
				while(getline(linestream, aLinePart, '\t')){
					lineParts.push_back(aLinePart);
				}
				if(lineParts.size() >= 11){
					// readID == [0], refID == [2], start == [3], cigar == [5], readSeq == [9]
					
					// Ignore reads with Ns
					if(lineParts[9].find("N") == string::npos){
						
						// Process cigar string
						stringstream startstream(lineParts[3]);
						int alignStart;
						startstream >> alignStart;
						alignStart = alignStart - 1;
						string fixedSeq;
						int valStartI = 0;
						int seqPointI = 0;
						for(int i=0; i<lineParts[5].length(); i++){
							if(!isdigit(lineParts[5].at(i))){
								
								char action = lineParts[5].at(i);
								string valuestr = lineParts[5].substr(valStartI, i-valStartI);
								stringstream valuess(valuestr);
								int value;
								valuess >> value;
								
								switch (action){
									case('M'):  // Match
										fixedSeq.append(lineParts[9].substr(seqPointI,value));
										seqPointI += value;
										break;
										
									case('S'):  // Soft-trim
										if(seqPointI == 0){
											//alignStart += value;
											seqPointI += value;
										}
										break;
										
									case('N'): { // Intron
										int alignEnd = alignStart + fixedSeq.length() - 1;
										reads.push_back(AlignedRead(fixedSeq, alignStart, alignEnd));
										alignStart = alignEnd + 1 + value;
										fixedSeq.clear();
										count++;
										}
										break;
										
									case('D'):  // Del
										for(int i=0; i<value; i++){
											fixedSeq.push_back(' ');
										}
										break;
										
									case('I'):  // Ins
										seqPointI += value;
										break;
								}
								
								valStartI = i+1;
							}
						}
						if(!(fixedSeq.empty())){
							int alignEnd = alignStart + fixedSeq.length() - 1;
							reads.push_back(AlignedRead(fixedSeq, alignStart, alignEnd));
							count++;
						}
					}
				}
			}
		}
	}
	catch(const boost::iostreams::gzip_error& e) {
		cerr << "Error while reading SAM file " << inSamFileName << endl;
		cerr << e.what() << endl;
	}
	fileifs.close();
	infile.reset();
	
	std::sort(reads.begin(), reads.end());
	
	return count;
}

	/*** Test reads from all samples over a SNP coord and print results if suitable **/
bool SNPTallyer::testSNP(const unsigned int snpCoord, const char refBase, const string& refID){
	bool printed = false;
	
	/* Base tally per sample i: 0 = A, 1 = T, 2 = C, 3 = G */
	int refBaseI = -1;
	switch(refBase){
		case 'A':
		case 'a':
			refBaseI = 0;
			break;
		case 'T':
		case 't':
			refBaseI = 1;
			break;
		case 'C':
		case 'c':
			refBaseI = 2;
			break;
		case 'G':
		case 'g':
			refBaseI = 3;
	}
	
	vector< vector<unsigned int> > baseTally;
	vector<unsigned int> totalReads (numSamples, 0);
	for(int sNum=0; sNum < numSamples; sNum++){
		baseTally.push_back(vector<unsigned int>(4, 0));
	}
	
	#pragma omp parallel for
	for(int sNum=0; sNum < numSamples; sNum++){	
		totalReads[sNum] = tallyBases(snpCoord, reads[sNum], baseTally[sNum]);
		//cout << totalReads[sNum] << " reads over snp " << refID << " " << snpCoord << " for sample " << sNum << endl;
	}
	
	vector<bool> basesToPrint(4, false);
	for(int baseI=0; baseI < 4; baseI++){
		if(baseI != refBaseI){
			bool printBase = false;
			for(int sNum=0; sNum < numSamples && !printBase; sNum++){
				if(baseTally[sNum][baseI] >= readDepthMin){
					printBase = true;
				}
			}
			if(printBase){
				printed = true;
				basesToPrint[baseI] = true;
			}
		}
	}
	if(printed){
		switch(outFormat){
			case 1:
				for(int baseI=0; baseI < 4; baseI++){
					if(basesToPrint[baseI]){	
						outtabfile << refID << "\t" << snpCoord << "\t" << refBase << "\t";
						switch(baseI){
							case 0:
								outtabfile << 'A';
								break;
							case 1:
								outtabfile << 'T';
								break;
							case 2:
								outtabfile << 'C';
								break;
							case 3:
								outtabfile << 'G';
						}
						for(int sNum=0; sNum < numSamples; sNum++){
							unsigned int numOther = totalReads[sNum] - baseTally[sNum][baseI];
							outtabfile << "\t" << baseTally[sNum][baseI] << "\t" << numOther;
						}
						outtabfile << "\n";
					}
				}
				break;
			case 2:
				for(int baseI=0; baseI < 4; baseI++){
					if(basesToPrint[baseI] || baseI == refBaseI){
						outtabfile << refID << "\t" << snpCoord << "\t";
						switch(baseI){
							case 0:
								outtabfile << 'A';
								break;
							case 1:
								outtabfile << 'T';
								break;
							case 2:
								outtabfile << 'C';
								break;
							case 3:
								outtabfile << 'G';
						}
						if(baseI == refBaseI){
							outtabfile << '*';
						}
						for(int sNum=0; sNum < numSamples; sNum++){
							outtabfile << "\t" << baseTally[sNum][baseI];
						}
						outtabfile << "\n";
					}
				}
				break;
			case 3:
				outtabfile << refID << "\t" << snpCoord << "\t" << refBase;
				for(int sNum=0; sNum < numSamples; sNum++){
					outtabfile << "\t" << baseTally[sNum][0];
					outtabfile << "\t" << baseTally[sNum][1];
					outtabfile << "\t" << baseTally[sNum][2];
					outtabfile << "\t" << baseTally[sNum][3];
				}
				outtabfile << "\n";
				break;
		}
		
	}
	return printed;
}

unsigned int SNPTallyer::tallyBases(const unsigned int snpCoord, const vector<AlignedRead>& sampReads, vector<unsigned int>& baseTally){
	/* Base tally : 0 = A, 1 = T, 2 = C, 3 = G */
	unsigned int totalReads = 0;
	
	// Narrow down search range
	int searchStart = 0;
	int searchEnd = sampReads.size()-1;
	bool change = true;
	while(change && searchStart < searchEnd){		
		change = false;
		int searchMid = searchStart + ((searchEnd - searchStart) / 2);
		if(sampReads[searchMid].start() > snpCoord){
			change = true;
			searchEnd = searchMid - 1;
		}else if(sampReads[searchMid].end()+readEndBuffer < snpCoord){
			change = true;
			searchStart = searchMid + 1;
		}
	}
	
	// Search of reads over SNP coord
	for(int i=searchStart; i<=searchEnd; i++){
		if(sampReads[i].start() > snpCoord){
			break;
		}else if(sampReads[i].end() >= snpCoord){
			if(sampReads[i].start()+edgeBuffer <= snpCoord && sampReads[i].end()-edgeBuffer >= snpCoord){
				totalReads++;
				const unsigned int readSNPCoord = snpCoord - sampReads[i].start();
				switch(sampReads[i][readSNPCoord]){
					case 'A':
					case 'a':
						baseTally[0] += 1;
						break;
					case 'T':
					case 't':
						baseTally[1] += 1;
						break;
					case 'C':
					case 'c':
						baseTally[2] += 1;
						break;
					case 'G':
					case 'g':
						baseTally[3] += 1;
				}
			}
		}
	}
	
	return totalReads;
}
