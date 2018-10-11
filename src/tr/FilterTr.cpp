/*
 * FilterTr.cpp
 * 
 *  Created on: Dec 14, 2012
 *      Author: Hani Zakaria Girgis, PhD
 */

// Delete start
#include <iostream>
#include <string>
#include <algorithm>
// Delete end

#include "FilterTr.h"
#include "TrKVisitor.h"
#include "TrCsVisitor.h"
#include "TrSineVisitor.h"
#include "../utility/Util.h"
#include "../utility/TSD.h"
#include "../utility/EmptyTSD.h"
#include "../utility/TailFinder.h"
#include "../utility/Tail.h"
#include "../utility/EmptyTail.h"
#include "../utility/GlobAlignE.h"
#include "../utility/LocAlign.h"
#include "../exception/InvalidStateException.h"
#include "../exception/InvalidInputException.h"

using namespace std;
using namespace utility;
using namespace exception;

namespace tr {

FilterTr::FilterTr(string nameIn,const string* seqIn, vector<BackwardTr*>* bListIn, int kIn, string bedFileNameIn, int ltrIdIn, int minIn, int maxIn, int minLtrLenIn,int maxLtrLenIn) {

	seq = seqIn;
	cSeq = seq->c_str();
	bList = bListIn;
	name = nameIn;
	k = kIn;
	teList = new vector<LtrTe *>();
	canUseLtr = false;
	canUseLength = true;
	canUseSine = false; //changed
	canUseDNA = true;
	canUsePpt = false;
	canUseTsd = false;//changed
	bedFileName = bedFileNameIn;

	init = (int) 'N';
	ltrId = ltrIdIn;
   
	tsdW = 100;
	tsdT = 4;

	tailW = 500;
	tailT = 12; //10
	min = minIn;
	max = maxIn;
	minLtrLen = minLtrLenIn;
	maxLtrLen = maxLtrLenIn;
    tightenBounds();
	//adjust();
	filter();
}


FilterTr::~FilterTr() {
	Util::deleteInVector(teList);	
	teList->clear();
	delete teList;
}

void FilterTr::tightenBounds()
{   //cout<<"Tightening bounds"<<endl;

	int chromLen = seq->length();

	vector<BackwardTr*> * newList = new vector<BackwardTr*>();

	for (int i=0;i<bList->size();i++)
	{
		BackwardTr * ltr = bList->at(i);

		int s1 = ltr->getS1();
		int e1 = ltr->getE1();
		int upLen = e1+k-s1;

		int s2 = ltr->getS2();
		int e2 = ltr->getE2();
		int downLen = e2+k-s2;

		int midpoint = e1+((s2-e1)/2);
		//cout <<"midpoint= "<<midpoint<<endl;

		//cout << "Before adjustment" << endl;
		//cout << "S1->" << s1 << " E1->" << e1 << " S2->" << s2 << " E2->" << e2 << endl;

		int upDiff = maxLtrLen - upLen;
		int downDiff = maxLtrLen - downLen;
		
		
		upDiff = upDiff < 0 ? 0 : upDiff;
		downDiff = downDiff < 0 ? 0: downDiff;

		int upStart = s1 - upDiff;
		int upEnd = upStart + upLen + 2 *upDiff + k;

		int downStart = s2 - downDiff;
		int downEnd = downStart + downLen + 2*downDiff + k;

		if (upStart >=0 && upEnd <chromLen && downStart >= 0 && downEnd<chromLen){

			//cout << "updiff=" << upDiff << endl;
			//cout << "downDiff=" << downDiff << endl;
			//cout << "upslice=" << upStart << ":" << upEnd << endl;

			upLen = upEnd > midpoint ? midpoint - upStart : upEnd-upStart; // Does alignment window encroach on midpoint between LTRs?

			//cout << "upslice=" << upStart << ":" << upStart+upLen<< endl;

			string upstream = seq->substr(upStart, upLen);
			//cout<<convertNucleotides(upstream)<<endl;

			downStart = downStart < midpoint ? midpoint : downStart;

			//cout << "downslice=" << downStart << ":" << downEnd<< endl;

			string downstream = seq->substr(downStart, downEnd-downStart);
			//cout<<convertNucleotides(downstream)<<endl;

			LocAlign *local = new LocAlign(upstream.c_str(), 0, upstream.length(), downstream.c_str(), 0, downstream.length(), 2, -3, 5, 2);
			//  cout<<"top "<<local->getQueryStart()<<":"<<local->getQueryEnd()<<endl;
			int newS1 = upStart + local->getQueryStart();
			int newE1 = upStart + local->getQueryEnd();
			//cout<<"bottom "<<local->getReferenceStart()<<":"<<local->getReferenceEnd()<<endl;
			int newS2 = downStart + local->getReferenceStart();
			int newE2 = downStart + local->getReferenceEnd();

			//cout<<"After adjustment"<<endl;
			//cout<<"S1->"<<newS1<<"E1->"<<newE1<<" S2->"<<newS2<<" E2->"<<newE2<<endl;
			double id = local->getIdentity() * 100;
			//cout<<"ID="<<id<<endl;
			BackwardTr *updated = new BackwardTr(newS1, newE1, newS2, newE2);

			if (id >= ltrId)
			{
				newList->push_back(updated);
			}

			//bList->at(i) = new BackwardTr(newS1,newE1,newS2,newE2);

			bool overlaps = true;

			int j = i + 1;

			// Avoiding duplicates
			while (overlaps && j < bList->size())
			{

				BackwardTr *next = bList->at(j);

				int nextStart = next->getS1();

				if (nextStart > updated->getE1())
				{
					//cout<<"no overlap"<<endl;
					overlaps = false;
					i = j - 1;
				}
				else
				{
					//cout<<"overlaps"<<endl;
					j += 1;
				}
			}
		}
		

	}

	//Removing duplicates
    std::sort( newList->begin(),newList->end(),[](BackwardTr * first, BackwardTr * second){ return first->getS1() < second->getS1();});
	//cout << "BLIST_SIZE === " << bList->size() << endl;

	Util::deleteInVector(bList);
	//bList->clear();

	
	for( int i = 0; i<(int)newList->size()-1;i++){

		BackwardTr * first = newList->at(i);
		int j = i+1;

		BackwardTr * second = newList->at(j);


			while (second->getS1() < first->getE1() && j<newList->size()-1)
		{

			j+=1;
			second = newList->at(j);

		}
		
		//cout<<"pushing back first"<<first->toString()<<endl;
		bList->push_back(first);
		//cout<<"pushing back second"<<second->toString()<<endl;
		bList->push_back(second);
		i = j;


	}

	//cout<<"BLIST_SIZE === "<<bList->size()<<endl;
	//bList = newList;

	//cout << "BLIST_SIZE === " << bList->size() << endl;

	delete newList;
	
}



void FilterTr::adjust() {
	int seqEnd = seq->size() - 1;
	TrKVisitor * kVisitor = new TrKVisitor(k, seqEnd); 
	int size = bList->size();
	for (int i = 0; i < size; i++) {
		//cout << bList->at(i)->toString()<<endl;
		bList->at(i)->accept(kVisitor); //adds k to index of each end TR
		//cout << bList->at(i)->toString()<<endl;
	}
}


void FilterTr::filter() {
	if (canUseLtr) {
		//cout<<"Entering alignment filter"<<endl;
		filterAcc2Ltr();
	} else {
		fillTeList();
		//cout<<"I shouldn't see this"<<endl;
	}
   // cerr<<"Past LTR"<<endl;
	if (canUseSine) {
		filterAcc2Sine();
	}
	
	if(canUseLength)
	{
	
		filterAcc2Length();
	}
	//cerr << "Past Length" << endl;

	if (canUseDNA)
	{
		filterAcc2DNA();
	}
//	cerr << "Past DNA" << endl;

	if (canUsePpt)
	{
		filterAcc2Ppt();
	}
    //cerr<<"Past Ppt"<<endl;
	if (canUseTsd) {
		filterAcc2Tsd();
	}
	//cerr<<"Past Tsd"<<endl;
}

void FilterTr::fillTeList() {

	//cout<<"inside fillTeList"<<endl;
	if (teList->size() != 0)
	{ 
			string msg("The TE list must be empty. The current size is: ");
		msg.append(Util::int2string(teList->size()));
		msg.append(".");
		throw InvalidStateException(msg);
	}

	int size = bList->size();
	//cout<<"SIZE IS :"<<size<<endl;
	for (int i = 0; i < size; i++) {
		teList->push_back(
				new LtrTe(bList->at(i), EmptyTSD::getInstance(),
						EmptyTail::getInstance()));
	}//cout<<"leaving fillTeList"<<endl;
}

void FilterTr::filterAcc2Ltr() {
	if (teList->size() != 0) {
		string msg("The TE list mBackwardTr * ltr = bList->at(i);ust be empty. The current size is: ");
		msg.append(Util::int2string(teList->size()));
		msg.append(".");
		throw InvalidStateException(msg);
	}
	int size = bList->size();

	for (int i = 0; i < size; i++) {
		BackwardTr * ltr = bList->at(i);

		// Overlap between the two LTRs
		int s1 = ltr->getS1();
		int e1 = ltr->getE1();
		int s2 = ltr->getS2();
		int e2 = ltr->getE2();
		//cout<<"S1="<<s1<<" e2="<<e2<<endl;
		

		if (!Util::isOverlapping(s1, e1, s2, e2)) {

			
			// Distance between the two LTRs
			int sep = abs(s2 - e2 + 1);
			if (sep >= ltrSep) {
				// Length and identity between the two LTRs
				TrCsVisitor * v = new TrCsVisitor(cSeq, minLtrLen, ltrId);
				ltr->accept(v);

				if (v->getIsGood()) {
					//cout <<"kept"<<endl;
					//cout << ltr->toString() << endl;
					teList->push_back(
							new LtrTe(ltr, EmptyTSD::getInstance(),
									EmptyTail::getInstance()));

				}
				else{
					//cout<<"Not kept:"<<endl;
					//cout << ltr->toString() << endl;
				}
				
				delete v;
			}
		}
		
	}
	//cout << "TElist size=" << teList->size() << endl;
}

void FilterTr::filterAcc2Length(){

	vector<LtrTe *> *temp = new vector<LtrTe *>();
	int size = teList->size();

	for (int i = 0; i < size; i++)
	{
		LtrTe *te = teList->at(i);

		BackwardTr *ltr = te->getLtr();
		
		int upLen = ltr->getE1()- ltr->getS1();
		int downLen = ltr->getE2() - ltr->getS2();

		int total = ltr->getE2()-ltr->getS1();
		
		bool upFits = upLen >= minLtrLen && upLen <=2* maxLtrLen;
		bool downFits = downLen >= minLtrLen && downLen <=2* maxLtrLen;
		bool totalFits = total >= min && total <= max;

		if(upFits && downFits && totalFits){
			//cout << "Keeping LTR " << ltr->getS1() << " "<< ltr->getE2() << " length=" << total << endl;
			LtrTe * longEnough = new LtrTe(*te);
			temp->push_back(longEnough);
		}
		else{
			//cout<<"Deleting LTR "<<ltr->getS1()<<" "<<ltr->getE2()<<" length="<<total<<endl;
		
	}


}
	Util::deleteInVector(teList);
	teList->clear();
	teList = temp;
	
	
}
// gets rid of 
void FilterTr::filterAcc2Sine() {
	vector<LtrTe *> * temp = new vector<LtrTe *>();
	int size = teList->size();

	for (int i = 0; i < size; i++) {
		LtrTe * te = teList->at(i);
		
		BackwardTr * ltr = te->getLtr();
        //cout<<ltr->toString()<<endl;
		//cout << "==========================================================================" << endl;
		TrSineVisitor * sineV = new TrSineVisitor(seq, tailT, 100, tsdT);
		ltr->accept(sineV);
		bool isTwoSines = sineV->isTwoSines();
		delete sineV;

		if (isTwoSines == false) {
			LtrTe * teWoSine = new LtrTe(*te);
			temp->push_back(teWoSine);
		}
	}

	Util::deleteInVector(teList);
	teList->clear();
	teList = temp;
}
// Filters according to target side duplication
void FilterTr::filterAcc2Tsd() {
	vector<LtrTe *> * temp = new vector<LtrTe *>();
	int size = teList->size();

	for (int i = 0; i < size; i++) {
		LtrTe * te = teList->at(i);
		
		BackwardTr * ltr = te->getLtr();
		TSD * tsd = new TSD(seq, ltr, tsdW, init);
		if (tsd->getTsdSize() > tsdT) {
			LtrTe * teWithTsd = new LtrTe(*te);
			teWithTsd->setTsd(tsd);
			temp->push_back(teWithTsd);
		}
		delete tsd;

	}

	Util::deleteInVector(teList);
	teList->clear();
	teList = temp;
}
// Filter according to polypurine tail. Existence of ppt proves it is ltr
void FilterTr::filterAcc2Ppt() {

	vector<LtrTe *> * temp = new vector<LtrTe *>();
	int size = teList->size();

	for (int i = 0; i < size; i++) {
		LtrTe * te = teList->at(i);
		BackwardTr * ltr = te->getLtr();
		Location *location = new Location(ltr->getS1(), ltr->getE2());

		TailFinder * finder = new TailFinder(seq, location, TailFinder::MARK_P,
				tailW, tailT);		
		vector<int> * tailVec = finder->getTail();
		if (tailVec->size() == 4) {
			string strand;
			int strandInt = tailVec->at(3);
			if (strandInt == 1) {
				strand = "+";
			} else if (strandInt == -1) {
				strand = "-";
			} else {
				string msg("Invalid strand. ");
				msg.append("Valid strands are 1 and -1, but received: ");
				msg.append(Util::int2string(strandInt));
				msg.append(".");
				cerr <<"About to throw!!!!"<<endl;
				throw InvalidInputException(msg);
			}
			Tail * tail = new Tail(tailVec->at(0), tailVec->at(1), strand,
					tailVec->at(2));
			LtrTe * teWithTail = new LtrTe(*te);
			teWithTail->setPpt(tail);
			delete tail;
			temp->push_back(teWithTail);
		}

			delete finder;
		delete location;
	}

	Util::deleteInVector(teList);
	teList->clear();
	teList = temp;
}

void FilterTr::filterAcc2DNA(){

	//cout<<"Inside DNA filter"<<endl;

	vector<LtrTe *> *filtered = new vector<LtrTe *>();

	int size = teList->size();

	for (int i = 0; i < size; i++)
	{

		LtrTe *te = teList->at(i);
		BackwardTr *ltr = te->getLtr();
		//cout<<"Processing"<<ltr->toString()<<endl;

		//Is there a TIR in the left detection

		int s1 = ltr->getS1();
		int e1 = ltr->getE1();

		int len1 = e1-s1;

		//cout<<"OG size="<<len1<<endl;

		len1 = len1 < minLtrLen ? len1/2 : minLtrLen/2;

		string deleteThis = seq->substr(s1,len1);

		//out<<convertNucleotides(deleteThis)<<endl;

		const char * leftFirst = deleteThis.c_str();

		//cout<<"leftFirst size ="<<seq->substr(s1,len1).length()<<endl;

		//cout<<"one"<<endl;
      
		string * temp = new string();
		
		
		Util::revCompDig(cSeq,e1-len1,e1,temp);

		//cout << "rightFirst size="<<temp->length() << endl;

		//cout<<"two"<<endl;

		//cout <<convertNucleotides(*temp)<<endl;

		const char * rightFirstRC = temp->c_str();
		//cout<<"three"<<endl;

       // cout<<"leftFirst length ="<<deleteThis.length()-1<<endl;
		//cout <<"rightFirstRC length ="<<temp->length()-1<<endl;

		 LocAlign *firstAlign = new LocAlign(leftFirst, 0, deleteThis.length()- 1, rightFirstRC, 0, temp->length() - 1, 2, -3, 5, 2);
		//cout<<"3.5"<<endl;

		bool leftTIR = firstAlign->getLength() >= 15;
		//cout<<"Alignment length ="<<firstAlign->getLength()<<endl;
		
		//cout << "four" << endl;

		//Is there a TIR in the right detection

		int s2 = ltr->getS2();
		int e2 = ltr->getE2();

		int len2 = e2-s2;


		len2 = len2 < minLtrLen ? len2 / 2 : minLtrLen / 2;

		string deleteThis2 = seq->substr(s2,len2);

		//cout<<convertNucleotides(deleteThis2)<<endl;

		const char *leftSecond = deleteThis2.c_str();

		//cout << "five" << endl;

		string * temp2 = new string();

		Util::revCompDig(cSeq, e2 - len2, e2, temp2);

		//cout << "six" << endl;

		//cout<<convertNucleotides(*temp2)<<endl;

		const char *rightSecondRC = temp2->c_str();

		LocAlign * secondAlign = new LocAlign(leftSecond, 0, deleteThis2.length()-1, rightSecondRC, 0, temp2->length()-1, 2, -3, 5, 2);

		//cout << "seven" << endl;

		bool rightTIR = secondAlign->getLength() >= 15;
		//cout <<secondAlign->getLength()<<endl;

		bool twoDNATransposons = leftTIR && rightTIR;

        if(!twoDNATransposons){

			LtrTe * notDNA = new LtrTe(ltr,EmptyTSD::getInstance(),EmptyTail::getInstance());

				filtered->push_back(notDNA);
		}
		else{

			//cout<<ltr->toString()<<"  consists of two adjacent DNA transposons"<<endl;
		}


	}

		Util::deleteInVector(teList);
		teList->clear();
		teList = filtered;
}

vector<LtrTe *> * FilterTr::getTeList() {
	return teList;
}

void FilterTr::bedFormat(int start, int end)
{
	cout << bedFileName << endl;
	ofstream output;
	output.open(bedFileName);
	int size = teList->size();

	//cout<<"FINAL TE SIZE ="<<size<<endl;

	for (int i = 0; i < size; i++)
	{
		LtrTe *curr = teList->at(i);
		BackwardTr* ltr = curr->getLtr();

		//cout<<ltr->toString()<<endl;

		//cout << name << "\t" << ltr->getS1() << '\t' << ltr->getE2() << "\t" << ltr->getS1() << '\t' << ltr->getE1() << "\t" << ltr->getS2() << '\t' << ltr->getE2() << endl;

		output << name
			   << "\t" << ltr->getS1() << '\t' << ltr->getE2() << "\t" << ltr->getS1() << '\t' << ltr->getE1()<<"\t"<<ltr->getS2()<< '\t' << ltr->getE2() << endl;
		//output<<convertNucleotides(upstream)<<endl;
		//output<<convertNucleotides(downstream)<<endl;

		
	}
	output.close();
}

string FilterTr::convertNucleotides(string str){
    string ans ="";
	for (int i = 0;i<str.length();i++){
        
		int curr = (int)str.at(i);
        int bp;

		switch(curr){

			case 0:
			    bp = 'A';
				break;
			case 1:
			    bp = 'C';
				break;
			case 2:
				bp = 'G';
				break;
			case 3:
				bp = 'T';
				break;
			default:
			    bp = 'N';

		}

		ans+=bp;

	}

    return ans;

}
void FilterTr::fullFormat(int start, int end)
{   
	cout << bedFileName << endl;
	ofstream output;
	output.open(bedFileName);

	int step =20;

		for (int i = 0; i < teList->size(); i++)
	{
		LtrTe *curr = teList->at(i);
		BackwardTr *ltr = curr->getLtr();
		
		int s1 = ltr->getS1();
		int e1 = ltr->getE1();
		int len1 = e1-s1;
		int x = maxLtrLen-len1;

		int s2 = ltr->getS2();
		int e2 = ltr->getE2();
		int len2 = e2-s2;
		int y = maxLtrLen-len2;

		int bound = x<y ? x:y;

	    if(len2<minLtrLen && bound<0.75*maxLtrLen){

		GlobAlignE *align = new GlobAlignE(cSeq, s1, e1, cSeq, s2, e2, 1, -1, 4, 1);

		int id = 100* align->getIdentity();

        //expand downstream
		for (int i =0;i<bound;i+=step){
           // cout<<"expanding downstream!"<<endl;
            GlobAlignE * extension = new GlobAlignE(cSeq,e1,e1+step-1,cSeq,e2,e2+step-1,1,-1,5,2);
			int extensionId = extension->getIdentity() *100;

			if(abs(id-extensionId)>=10){
				e2+=step-1;
				e1+=step-1;
				cout << "chrX" << "	" << s1 << " " << e2 << endl;
				break;
			}

			else if (abs(id - extensionId) >= 10)
			{
				cout << "chrX" << "	" << s1 << " " << e2 << endl;
				break;
			}
		}
        //expand upstream
		for (int i = 0; i < bound; i += step)
		{
            cout<<"expanding upstream"<<endl;
			GlobAlignE *extension = new GlobAlignE(cSeq, s1-step, s1, cSeq, s2-step, s2, 2, -3, 4, 1);
			int extensionId = extension->getIdentity() * 100;

			if (abs(id - extensionId) >= 10 && extensionId>50)
			{
				s2 -=  step - 1;
				s1 -=  step - 1;
				cout << "chrX"<< "	" << s1 << " " << e2 << endl;
				break;
			}

			else if (abs(id - extensionId) >= 10)
			{
				cout << "chrX"<< " " << s1 << " " << e2 << endl;
				break;
			}
			
		}

		
		}

		if (e2 - s1 > minLtrLen)
		{
			output << "chrX" << "\t" << s1 << "\t" << e2<< endl;
		}
	}
	output.close();
}
} /* namespace tr */
