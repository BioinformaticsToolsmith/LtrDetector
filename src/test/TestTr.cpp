/*
 * TestTr.cpp
 *
 *  Created on: Nov 30, 2012
 *      Author: Hani Zakaria Girgis, PhD
 */

#include "../nonltr/ChromosomeOneDigit.h"
#include "../nonltr/Chromosome.h"
#include "../nonltr/ChromosomeRandom.h"

#include "../tr/ScorerTr.h"
#include "../tr/FilterTr.h"
#include "../tr/Tr.h"
#include "../tr/ForwardTr.h"
#include "../tr/BackwardTr.h"
#include "../tr/TrCollector.h"
#include "../utility/Util.h"
#include "../utility/LCS.h"
#include "../utility/LCSLen.h"
#include "../utility/Location.h"
#include "../utility/ILocation.h"
#include "../utility/LCSubStr.h"
#include "../utility/TSD.h"

#include <cstdlib>
#include <iostream>
#include <fstream>

using namespace std;

using namespace nonltr;
using namespace utility;
using namespace tr;



int main(int argc, char * argv[]) {
	
	//string chromFile = string("../src/FlyBaseR5.48/chrX-chromosome.fasta");
    //string csvFile = string("../src/output/trial1.csv");
	//string bedFile = string("../src/output/trial1.bed");
	//string chromFile = string(argv[1]);
	string chromDir = string(argv[1]);
	string outputDir = string(argv[2]);
	
	int minDist = atoi(argv[3]);
    int maxDist = atoi(argv[4]);

	int minLenLTR = atoi(argv[5]);
	int maxLenLTR = atoi(argv[6]);

	int identity = atoi(argv[7]);

	int k = atoi(argv[8]);
	int minPlateauLength = atoi(argv[9]);
	int gapTol =atoi(argv[10]);

	

	
	//string csvFile = string(argv[9]);
	
	//cout<<"FIRST ARGUMENT"<<chromFile<<endl;
	//cout<<"MIN ARGUMENT"<<minDist<<endl;
	//cout<<"MAX ARGUMENT"<<maxDist<<endl;
	//cout<<"LAST ARGUMENT"<<bedFile<<endl;;
    //exit(1);

	vector <string> * chromList = new vector<string>();
	cout<<"Reading chromosome directory"<<endl;
	Util::readChromList(chromDir, chromList,"fa");

	ChromosomeOneDigit *chrom;

	TrCollector *collector;
	
	#pragma omp parallel for schedule(dynamic) num_threads(3)
	
	 for (int i = 0; i < chromList->size(); i++)
	{

		string chromFile = chromList->at(i);
		cout<<"Processing: "<<chromFile<<endl;

		int nameBegin = chromFile.find_last_of("/")+1;
		int nameEnd = chromFile.find_last_of(".") ;
		int len = nameEnd - nameBegin ;

		string name = chromFile.substr(nameBegin,len);

		
		string bedFile = outputDir+"/"+name+"Detector.bed";

		//chrom = new ChromosomeOneDigit(chromFile);
		collector = new TrCollector(chromFile, bedFile,name, minDist, maxDist, minLenLTR, maxLenLTR,identity,k, minPlateauLength, gapTol);

		cout << "Output from" << name << " found in: " << bedFile << endl;

		delete collector;
		
	}

   
	//TrCollector * collector = new TrCollector(chrom, 14, 2000, 11000, 11000, 50, 10, 10, "../src/output/test1.csv", indexFile);
	//TrCollector *collector = new TrCollector(chrom, 14, 0, 0, 0);
	// TrCollector * collector = new TrCollector(chrom, 9, 2000, 11000, 1000);
	//vector<LtrTe *> * teList = collector->getTeList();
	//cout << "Total number of LTR TE is: " << teList->size() << endl;
	//collector->printIndex(indexFile);
	

	return 0;

}

