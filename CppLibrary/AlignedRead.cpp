#include <string>
#include "AlignedRead.h"
using namespace std;

/* By Andrew Spriggs, CSIRO Ag&Food, 2018 */
/* andrew.spriggs@csiro.au */
/* https://github.com/spriggsy83 */

AlignedRead::AlignedRead(const string& newSeq, int newStart, int newEnd){
	sequence = newSeq;
	alignedStart = newStart;
	alignedEnd = newEnd;
}

AlignedRead::AlignedRead(const AlignedRead& copySource){
	alignedStart = copySource.alignedStart;
	alignedEnd = copySource.alignedEnd;
	sequence = copySource.sequence;
}

AlignedRead& AlignedRead::operator= (const AlignedRead& copySource){
	if(this != &copySource){
		alignedStart = copySource.alignedStart;
		alignedEnd = copySource.alignedEnd;
		sequence = copySource.sequence;
	}
	return *this;
}

int AlignedRead::start() const{
	return alignedStart;
}

int AlignedRead::end() const{
	return alignedEnd;
}

const char& AlignedRead::operator[] (int i) const{
	return sequence[i];
}

string AlignedRead::getSeq() const{
	return sequence;
}

bool AlignedRead::operator < (const AlignedRead& otherRead) const{
	return (alignedStart < otherRead.alignedStart);
}

