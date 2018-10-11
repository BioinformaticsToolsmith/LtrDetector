/*
 * TrCollector.cpp
 *
 *  Created on: Jan 2, 2013
 *      Author: Hani Zakaria Girgis, PhD
 */

#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include "TrCollector.h"
#include "ScorerTr.h"
#include "DetectorTr.h"
#include "FilterTr.h"
#include "MatchTr.h"
#include "../nonltr/ChromosomeOneDigit.h"

#include "../utility/Util.h"

using namespace std;
using namespace utility;

namespace tr {

TrCollector::TrCollector(std::string chromFileIn, std::string bedFileNameIn, std::string nameIn, int minIn,
						 int maxIn, int ltrMinIn, int ltrMaxIn, int idIn, int kIn, int plateauLenIn, int gapTolIn)
{
	k = kIn;
	min = minIn;
	max = maxIn;
	//d = dIn;
	minPlateauLen = plateauLenIn;

	gapTol = gapTolIn;
	name = nameIn;
	ltrMin = ltrMinIn;
	ltrMax = ltrMaxIn;
	//csvFileName = csvFileNameIn;MatchTr
	bedFileName = bedFileNameIn;
    identity = idIn;
	chrom = new ChromosomeOneDigit(chromFileIn);

	teList = new vector<LtrTe *>();
	/*for (int y = min; y < max; y += d) {
		int e = y + d - 1;
		if (e > max) {
			e = max;
		}

		collect(y, e);
	}*/
	collect();

	// Sort detections according to the start site
	sort(teList->begin(), teList->end(), LtrTe::lessThan);

	// Testing start
	/*
	 int size = teList->size();
	 for (int i = 0; i < size; i++) {
	 cout << i << " " << teList->at(i)->toString() << endl;
chrX	445510	453014

	 }
	 */
	// Testing end
}

void TrCollector::collect() {

	//cerr << "Processing " << lower << ":" << upper << endl;
	//cerr << "Scoring ..." << endl;
    int maxScore = max - ltrMin <=2000 ? 2000: max-ltrMin;
	int minScore = min - ltrMax <= 2000 ? 2000 : min - ltrMax;
   // cerr <<"MINSCORE="<<minScore<<endl;
	ScorerTr * scorer = new ScorerTr(chrom, k, minScore,maxScore);
	int length = chrom->getEffectiveSize();
	
	//cout<<"Chrom length:"<<length<<endl;

	//cout<<"Passed scorerTR"<<endl;

	//cout << "passed matchtr" << endl;
	//scorer->scoresFormat(0,length);

	//cout<<"passed print scores"<<endl;
	scorer->scoresFormat(0,length);
	//scorer->bedFormat(0, chrom->size());
	int init = scorer->getInitialScore();
	//scorer->outputScores()

	MatchTr * matcher = new MatchTr(scorer->getScores(),k, init, bedFileName, min,max,ltrMin, minPlateauLen, gapTol,identity);

	//cout<<"passed matchtr"<<endl;

	//matcher->printFinalScores(0,length);

	//matcher->bedFormat(0,chrom->size());

	//cout<<"passed print bed"<<endl;

	//MatchTr * matcher = new MatchTr(scorer->getScores(), init, lower,upper, minPlateauLen,diffThres,gaptol);
	//DetectorTr * detector = new DetectorTr(scorer->getScores(), init);
    
	//cerr << "Filtering ..." << endl;
	FilterTr *filter = new FilterTr(name, chrom->getBase(), matcher->getRepeatCandidates(), k, bedFileName,identity,min,max,ltrMin,ltrMax);
	filter->bedFormat(0,chrom->size());
	
	//FilterTr * filter21427124 = new FilterTr(chrom->getSegment(),detector ->getBList(),k)
	//cerr << "Finished Filtering";
	//vector<LtrTe *> * part = filter->getTeList();
	//int size = part->size();

	/*for (int f = 0; f < size; f++) {
		teList->push_back(new LtrTe(*(part->at(f))));
	}

	delete filter;*/

	
	delete scorer;
}

TrCollector::~TrCollector() {


	Util::deleteInVector(teList);
	teList->clear();
	delete teList;
	delete chrom;
	
}

vector<LtrTe *> * TrCollector::getTeList() {
	return teList;
}

void TrCollector::printIndex(string outputFile) {
	ofstream outIndex;
	outIndex.open(outputFile.c_str(), ios::out /*| ios::app*/);

				  // Write the index of the repeat segment [x,y[ "exclusive" with respect with the start (chrK:start-end)
				  string header = chrom->getHeader();
	int size = teList->size();
	for (int j = 0; j < size; j++) {
		LtrTe * te = teList->at(j);

		// outIndex << header << ":" << te->getStart() << "-" << te->getEnd() /*+ 1*/
		// << endl;5000

		outIndex << te->toString(header) << endl;
	}
	outIndex.close();
}

void TrCollector::printMasked(string outputFile) {
	
	/*string baseCopy = string(*(chrom->getBase()));
	int size = teList->size();

	for (int j = 0; j < size; j++) {
		LtrTe * te = teList->at(j);
		int teStart = te->getStart();
		int teEnd = te->init, k, getEnd();

		for (int h = teStart; h <= teEnd; h++) {
			baseCopy[h] = tolower(baseCopy[h]);
		}
	}

	ofstream outMask;
	outMask.open(outputFile.c_str(), ios::out /*| ios::app);
	outMask << chrom->getHeader() << endl;
	int step = 50;
	int len = baseCopy.size();
	for (int i = 0; i < len; i = i + step) {
		int e = (i + step - 1 > len - 1) ? len - 1 : i + step - 1;
		for (int k = i; k <= e; k++) {
			outMask << baseCopy[k];
		}
		outMask << endl;
	}
	outMask.close();*/
}

} /* namespace tr */
