/*
 * Assembler.cpp
 *
 *  Created on: Nov 27, 2012
 *      Author: Hani Zakaria Girgis, PhD
 */

#include "TailFinder.h"

#include "../exception/InvalidInputException.h"
#include "../exception/InvalidStateException.h"
#include "../utility/Util.h"

#include <vector>

// #include <iostream>

using namespace std;
using namespace exception;

namespace utility {
TailFinder::TailFinder(const string * seqIn, ILocation * locIn, int whichTailIn,
		int winIn, int minLenIn) {
	seq = seqIn;
	loc = locIn;
	whichTail = whichTailIn;
	win = winIn;
	minLen = minLenIn;

	seedLen = 2;
	gapLen = 4;

	findMark();
}


TailFinder::~TailFinder() {
	tail->clear();
	delete tail;
}

void TailFinder::findMark() {

	int start = loc->getStart();
	
	int end = loc->getEnd();

	int len = loc->getLength();

	string * detection = new string(seq->begin() + start,
			seq->begin() + end + 1);

	int halfEnd = detection->size() / 2;

	halfEnd = (halfEnd > win) ? win : halfEnd;

	vector<int> * pstvTail = new vector<int>(); //stores info on positive tail

    cout<<"Repeat:"<<endl;
	// Test start
	for (int h = 0; h < detection->size(); h++)
	{
		cout << static_cast<int>(detection->at(h)) << " ";
	}
	cout << endl;
	// Test end

    cout<<"Tail Details:"<<endl;
	if (whichTail == MARK_A) {
		cout<<"Entering MarkA"<<endl;
		findMarkA(detection, pstvTail, 0, halfEnd);
	} else if (whichTail == MARK_P) {
		findMarkP(detection, pstvTail, 0, halfEnd);
	} else {
		string msg("Invalid mark. Valid marks are: ");
		msg += Util::int2string(MARK_A);
		msg += string(" and ");
		msg += Util::int2string(MARK_P);
		msg += string(". Received: ");
		msg += Util::int2string(whichTail);
		msg += ".";
		throw InvalidInputException(msg);
	}

	

			//(start,end,ratio)
			if (pstvTail->size() == 3)
			{
				(*pstvTail)[0] = start + pstvTail->at(0);
				(*pstvTail)[1] = start + pstvTail->at(1);
			}
			//string output = prettyFormatChrom(detection);
			//cout <<output <<endl;



			vector<int> *ngtvTail = new vector<int>(); //stores info on reverse complement tail
			string *rcDetection = new string();
			Util::revCompDig(detection, rcDetection);
			delete detection;

			if (whichTail == MARK_A)
			{   cout<<"MarkA for reverse complement"<<endl;
				findMarkA(rcDetection, ngtvTail, 0, halfEnd);
			}
			else if (whichTail == MARK_P)
			{
				findMarkP(rcDetection, ngtvTail, 0, halfEnd);
			}
			delete rcDetection;

			if (ngtvTail->size() == 3)
			{
				int ss = start + len - ngtvTail->at(0) - 1; //changed from at(1) to reflect alteration of BackwardTr
				cout << "ss=" << ss+1 << endl;
				int ee = start + len - ngtvTail->at(1) - 1; // changed from at(0)
				cout << "ee=" << ee << endl;
				(*ngtvTail)[0] = ss;
				(*ngtvTail)[1] = ee;
			}

			// The reverse sign is correct
			if (pstvTail->size() == 3 && ngtvTail->size() == 0)
			{
				tail = pstvTail;
				tail->push_back(-1);
				delete ngtvTail;
			}
			else if (pstvTail->size() == 0 && ngtvTail->size() == 3)
			{
				tail = ngtvTail;
				tail->push_back(1);
				delete pstvTail;
			}
			else if (pstvTail->size() == 0 && ngtvTail->size() == 0)
			{
				tail = pstvTail;
				delete ngtvTail;
			}
			else if (pstvTail->size() == 3 && ngtvTail->size() == 3)
			{ //if both exist, use the longer one
				int pstvLen = pstvTail->at(1) - pstvTail->at(0) + 1;
				int ngtvLen = ngtvTail->at(1) - ngtvTail->at(0) + 1;
				if (pstvLen > ngtvLen)
				{
					tail = pstvTail;
					tail->push_back(-1);
					delete ngtvTail;
				}
				else
				{
					tail = ngtvTail;
					tail->push_back(1);
					delete pstvTail;
				}
			}
			else
			{
				string msg = string("The tail must have three coordinates only. ");
				msg += string("The +ve tail has ");
				msg += Util::int2string(pstvTail->size());
				msg += string(" coordinates. The -ve tail has ");
				msg += Util::int2string(ngtvTail->size());
				msg += string(" coordinates.");
				throw InvalidStateException(msg);
			}

			// For testing		int winIn, int minLenIn) {

			/*
	 cout << start << "-" << end << endl;

	 if (tail->size() == 4) {
	 cout << ">> Tail length is: " << tail->at(1) - tail->at(0) + 1
	 << " start: " << tail->at(0) << " end: " << tail->at(1)
	 << " percentage is: " << tail->at(2) << "% strand is: "
	 << tail->at(3) << endl;
	 for (int i = start; i <= end; i++) {
	 if (i == tail->at(0)) {
	 cout << "_";
	 }

	 cout << (int) seq->at(i);
	 if (i == tail->at(1)) {
	 cout << "_";		int winIn, int minLenIn) {

	 }
	 }
	 cout << endl;

	 } else {
	 cout << "No PPT was found" << endl;
	 for (int i = start; i <= end; i++) {
	 cout << (int) seq->at(i);
	 }
	 cout << endl;
	 }
	 cout << endl;
	 */
			// End testing
}
string utility::TailFinder::prettyFormatChrom( string * detection){
    std:: string ans;
	for( int i =0;i<detection->length();i++){
		ans+= (int)detection->at(i);
	}
	return ans;
}

//used to find poly(A)
void TailFinder::findMarkA(string * detection, vector<int> * tail, int segStart,
		int halfEnd) {
	int code = 3;//code for T nucleotide, which is complement of A in mRNA
	double ratio = 0.8;

	// Find the first seed available then immediately break
	bool isSeedFound = false;                                                                
	int seedEnd;

	for (int y = 0; y <= halfEnd - seedLen + 1; y++) { 
		int count = 0;
		for (int x = y; x <= y + seedLen - 1; x++) {
			if (detection->at(x) == code) {
				count++;
			}
		}

		if (count == seedLen) {
			seedEnd = y + seedLen - 1;
			isSeedFound = true;
			cout<<"Seed found at"<<seedEnd-seedLen +1<<","<<seedEnd<<endl;
			break;
		}
	}
    
	// Extend initial seed, allowing for gaps of up to gapLen
	if (isSeedFound) {
		int gapCounter = 0;
		int i = seedEnd + 1;
		for (; i < halfEnd; i++) {
			if (static_cast<int>(detection->at(i)) == code) {
				gapCounter = 0;
			} else {
				gapCounter++;
			}

			if (gapCounter == gapLen) {
				break;
			}
		}
        
		int tailStart = seedEnd - seedLen + 1;
		int tailEnd = i;
		//cout<<"start:"<<tailStart<<" end:"<<tailEnd<<endl;
		//retreats back to last A nucleotide
		for (int j = i; j >= i - gapLen; j--) {
			if (static_cast<int>(detection->at(j)) == code) {
				tailEnd = j;
				break;
			}
		}
        //must end with T
		if (detection->at(tailEnd) != code) {
			string msg("Invalid tail. The tail does not end with T.");
			throw InvalidStateException(msg);
		}
        //cout<<"tailStart:"<<tailStart<<" tailEnd:"<<tailEnd<<endl;
		
		double pTailLen = tailEnd - tailStart + 1;
		double tCount = 0;
		for (int h = tailStart; h <= tailEnd; h++) {
			if (static_cast<int>(detection->at(h)) == code) {
				tCount++;
			}
		}
        /*cout<<" Longest A sequence is: "<<pTailLen<<endl;
		cout<<"tCount="<<tCount<<endl;
		cout<<"ratio ="<<(tCount/pTailLen)<<endl;*/

		if (pTailLen >= minLen && (tCount / pTailLen) >= ratio) {
			tail->push_back(tailStart + segStart); //tail[0] = start
			tail->push_back(tailEnd + segStart); //tail[1] = end
			double tCount = 0;
			cout<<"tCount1="<<tCount<<endl;
	

			tail->push_back(100 * tCount / pTailLen);
		}
/*
		if (tail->size() == 0 && halfEnd - segStart + 1 >= minLen) { //there is enough left to fit in a minLen tail
			int shift = segStart + seedEnd; // increment by 
			string * rest = new string(detection->begin() + seedEnd,
					detection->begin() + detection->size());
			findMarkA(rest, tail, shift, halfEnd - seedEnd);
			delete rest;
		}*/
	}
}

//used to find purines (A and G)
void TailFinder::findMarkP(string * detection, vector<int> * tail, int segStart,
		int halfEnd) {
	double ratio = 0.65;
	int codeC = 1;
	int codeT = 3;

	// Find a seed
	bool isSeedFound = false;
	int seedEnd;
	for (int y = 0; y <= halfEnd - seedLen + 1; y++) {
		int count = 0;
		for (int x = y; x <= y + seedLen - 1; x++) {
			if (detection->at(x) == codeC || detection->at(x) == codeT) {
				count++;
			}
		}

		if (count == seedLen) {
			seedEnd = y + seedLen - 1;
			isSeedFound = true;
			break;
		}
	}

	// Extend
	if (isSeedFound) {
		int gapCounter = 0;
		int i = seedEnd + 1;
		for (; i < halfEnd; i++) {
			if (detection->at(i) == codeC || detection->at(i) == codeT) {
				gapCounter = 0;
			} else {
				gapCounter++;
			}
			if (gapCounter == gapLen) {
				break;
			}
		}

		int tailStart = seedEnd - seedLen + 1;
		int tailEnd = i;
		for (int j = i; j >= i - gapLen; j--) {
			if (detection->at(j) == codeC || detection->at(j) == codeT) {
				tailEnd = j;
				break;
			}
		}

		if (!(detection->at(tailEnd) == codeC || detection->at(tailEnd) == codeT)) {
			string msg("The tail does not end with C or T. ");
			msg.append("The invalid base is: ");
			msg.append(Util::int2string((int) detection->at(tailEnd)));
			msg.append(".");
			throw InvalidStateException(msg);
		}

		double pTailLen = tailEnd - tailStart + 1;
		double ctCount = 0;
		for (int h = tailStart; h <= tailEnd; h++) {
			if (detection->at(h) == codeC || detection->at(h) == codeT) {
				ctCount++;
			}
		}

		if (pTailLen >= minLen && (ctCount / pTailLen) >= ratio) {
			tail->push_back(tailStart + segStart);
			tail->push_back(tailEnd + segStart);

			double pCount = 0;
			for (int g = tailStart; g <= tailEnd; g++) {
				if (detection->at(g) == codeC || detection->at(g) == codeT) {
					pCount++;
				}
			}
			tail->push_back(100 * pCount / pTailLen);
		}

		if (tail->size() == 0 && halfEnd - segStart + 1 >= minLen) {
			int shift = segStart + seedEnd;
			string * rest = new string(detection->begin() + seedEnd,
					detection->begin() + detection->size());
			findMarkP(rest, tail, shift, halfEnd - seedEnd);
			delete rest;
		}
	}
}

/*
 * The size of the tail indicates the following:
 * 0: no tail is found
 * 4: start, end, strand: 1 indicates pstv and -1 indicates ngtv.
 */
vector<int> * TailFinder::getTail() {
	return tail;
}

/**
 * If the vector representing the tail has a size of zero,
 * no tail is found.
 */
bool TailFinder::isTailFound() {
	bool r = tail->size() == 0 ? false : true;
	return r;
}

}
