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
#include <algorithm>
#include <vector>
#include <fstream>

using namespace std;

using namespace nonltr;
using namespace utility;
using namespace tr;

int getOptionalArg(vector<string> * args, string option){
    string answer = "-1";
	auto it = std::find(args->begin(),args->end(),option);

	if(it!=args->end() & ++it !=args->end()){

		answer = * it;
	}

	return stoi(answer);
}



int main(int argc, char * argv[]) {
	
	//string chromFile = string("../src/FlyBaseR5.48/chrX-chromosome.fasta");
    //string csvFile = string("../src/output/trial1.csv");
	//string bedFile = string("../src/output/trial1.bed");
	//string chromFile = string(argv[1]);

	int minDist = 2000;
	int maxDist = 18000;
	int minLenLTR = 100;
	int maxLenLTR = 2000;
	int identity = 85;
	int k = 14;
	int minPlateauLength = 10;
	int gapTol = 200;
	string chromDir ="";
	string outputDir = "";
	bool printRawScores = false;
	bool printCleanScores = false;

	//string options [] = {"-minLen", "-maxLen","-minLenLTR", "-maxLenLTR", "-id", "-k", "-plateauSeed","-gapTol"};

	std::vector<string> * argList = new vector<string>();

	for(int i = 1;i<argc;i++){

		argList->push_back(string(argv[i]));
	}
	
	auto src = std::find(argList->begin(),argList->end(),"-chromDir");

	if(src!=argList->end() & ++src !=argList->end()){

		 chromDir = *src;
	}

	else{
		cout<<" \'-chromDir\' is a required argument"<<endl;
		exit(1);
	}

	auto dest = std::find(argList->begin(),argList->end(),"-destDir");

	if(dest!=argList->end() & ++dest != argList->end()){

		outputDir = *dest;
	}
	else{
		cout<<" \'-destDir\' is a required argument"<<endl;
		exit(1);
	}

	int test = getOptionalArg(argList, "-minLen");

	if(test!=-1){
		minDist = test;
	}

	test = getOptionalArg(argList, "-maxLen");

	if(test!=-1){
		maxDist = test;
	}

	test = getOptionalArg(argList,"-minLenLTR");

	if(test!=-1){

		minLenLTR = test;
	}

	test = getOptionalArg(argList, "-maxLenLTR");

	if(test!=-1){
		maxLenLTR = test;
	}

	test = getOptionalArg(argList,"-id");

	if(test!= -1){

		identity = test;
	}

	test = getOptionalArg(argList, "-k");

	if(test!=-1){
		k = test;
	}

	test = getOptionalArg(argList, "-plateauSeed");

	if(test!=-1){
		minPlateauLength = test;
	}

	test = getOptionalArg(argList, "-gapTol");

	if(test!=-1){
		gapTol = test;
	}

	
	
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

		cout<<chromFile<<endl;
		cout<<bedFile<<endl;
		cout<<minDist<<endl;
		cout<<maxDist<<endl;
		cout<<minLenLTR<<endl;
		cout<<maxLenLTR<<endl;
		cout<<identity<<endl;
		cout<<k<<endl;
		cout<<minPlateauLength<<endl;
		cout<<gapTol<<endl;
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

