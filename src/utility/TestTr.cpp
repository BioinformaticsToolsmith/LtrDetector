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
#include "../tr/DetectorTr.h"
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

#include <iostream>
#include <fstream>

using namespace std;

using namespace nonltr;
using namespace utility;
using namespace tr;

void printScores(string outputFile, vector<int>* scores,
		ChromosomeOneDigit * chrom) {
	ofstream outScores;
	outScores.open(outputFile.c_str(), ios::out);

	int step = 50;

	outScores << chrom->getHeader() << endl;
	int len = scores->size();
	for (int i = 0; i < len; i = i + step) {
		int e = (i + step - 1 > len - 1) ? len - 1 : i + step - 1;
		for (int k = i; k <= e; k++) {
			outScores << scores->at(k) << " ";
		}
		outScores << endl;
	}
	outScores << endl;

	outScores.close();
}

int main(int argc, char * argv[]) {
	/*
	 string s1(
	 "323331033213003103021003332330002031221233231203103200303333001212002102303203330030032032232330200000033000003000003200000033002221202312230200001013");
	 string s2(
	 "033323323331221020103320000203230112010003300003020300300320310000030313333320003101002333322120231223010003230000311100030020020133301312332023333323");
	 LCS * lcs = new LCS(s1.c_str(), 0, 149, s2.c_str(), 0, 149);

	 cout << "The length of the LCS is: " << lcs->getLenCS() << endl;
	 lcs->printLcs();

	 // Read sequence file
	 cout << "Reading chromosom sequence" << endl;
	 */

	// CBB
	// string chromFile = string("/home/girgishz/myGenome/hg19/fa/chr20.fa");
	// string chromFile("/panfs/pan1/msdetector/hg19/chr20_upper_shuffled.fa");
	// string chromFile("/panfs/pan1/msdetector/repeatsDetector/hgTest/chr20_shuffled_ltr100.fa");
	// Zakarota
	// HG
	/*
	 string chromFile = string("/Users/zakarota/Data/Genomes/Hg19/chr20.fa");
	 string indexFile = string("/Users/zakarota/Data/HgTest/Ltr/chr20Ltr.coor");
	 */

	// MD
	string chromFile = string(
			"/Users/zakarota/Data/Repeats/FlyBaseR5.48/chrX-chromosome.fasta");

	// string indexFile = string("/Users/zakarota/Data/DmTest/Ltr/chrXLtr.coor");
	string indexFile = string("/Users/zakarota/Data/DmTest/Ltr/chrXTest.coor");

	const string randChromFile = string(
			"/Users/zakarota/Data/Repeats/FlyBaseR5.48/chrX-rand-chrom.fasta");

	// Test the generation of random chromosome
	/*
	 Chromosome * oChrom = new Chromosome(chromFile);
	 // Chromosome * oChrom = new Chromosome(randChromFile);
	 vector<char> * alpha = new vector<char>();
	 alpha->push_back('A');
	 alpha->push_back('C');
	 alpha->push_back('G');
	 alpha->push_back('T');
	 ChromosomeRandom * rChrom = new ChromosomeRandom(5, oChrom, 'N', alpha);

	 // Print the random chromosome
	 const string * base = rChrom->getBase();
	 const string header = rChrom->getHeader();
	 Util::writeFasta(*base, header, randChromFile);
	 */

	// Detect LTR
	ChromosomeOneDigit * chrom = new ChromosomeOneDigit(chromFile);
	TrCollector * collector = new TrCollector(chrom, 14, 2000, 11000, 11000);
	// TrCollector * collector = new TrCollector(chrom, 9, 2000, 3000, 1000);
	vector<LtrTe *> * teList = collector->getTeList();
	cout << "Total number of LTR TE is: " << teList->size() << endl;
	collector->printIndex(indexFile);
	delete collector;
	delete chrom;

	return 0;
}
