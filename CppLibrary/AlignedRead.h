#ifndef ALIGNEDREAD_H
#define ALIGNEDREAD_H

#include <string>
using namespace std;

/* By Andrew Spriggs, CSIRO Ag&Food, 2018 */
/* andrew.spriggs@csiro.au */
/* https://github.com/spriggsy83 */

class AlignedRead {
	string sequence;
	int alignedStart;
	int alignedEnd;
	
  public:
	AlignedRead(const string&, int, int);
	AlignedRead(const AlignedRead&);
	AlignedRead& operator= (const AlignedRead&);
	bool operator < (const AlignedRead& otherRead) const;
	int start() const;
	int end() const;
	const char& operator[] (int i) const;
	string getSeq() const;
};

#endif
