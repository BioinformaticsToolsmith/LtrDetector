/*
 * TrCollector.h
 *
 *  Created on: Jan 2, 2013
 *      Author: Hani Zakaria Girgis, PhD
 */
// 
#ifndef TRCOLLECTOR_H_
#define TRCOLLECTOR_H_

#include "LtrTe.h"
#include "../nonltr/ChromosomeOneDigit.h"

using namespace nonltr;

namespace tr {

class TrCollector {
private:
	ChromosomeOneDigit * chrom;
	int k;
	int min; // minimum separation distance between ltr
	int max; //absolute maximum ^
	int ltrMin;
	int ltrMax;
	//int d; // "delta" incremental separation between ltr on iteration
	int minPlateauLen;
	int diffThresh;
	
	int gapTol;
	int identity;
	std::string csvFileName;
	std::string bedFileName;
	std::string name;
	vector<LtrTe *> * teList;



	void collect();

public:
	TrCollector(std::string,std::string,std::string, int, int,int,int,int,int,int,int);
	virtual ~TrCollector();
	vector<LtrTe *> * getTeList();
	void printIndex(string);
	void printMasked(string);

};

} /* namespace tr */
#endif /* TRCOLLECTOR_H_ */
