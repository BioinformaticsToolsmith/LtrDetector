/*
  Created by Joseph V
  encia 21 February 2018
*/

#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>
#include <stack>

#include "MatchTr.h"
#include "ForwardTr.h"
#include "Candidate.h"
#include "PairContainer.h"
#include "../utility/Util.h"
#include "../exception/InvalidInputException.h"

using namespace std;

//using namespace nonltr;
using namespace utility;
using namespace exception;

namespace tr{

MatchTr::MatchTr(vector<int> * scoreListIn,int kIn, int initValueIn,string bedFileIn,int minIn, int maxIn, int ltrMinIn,int plateauLenIn,  int gapTolIn, int id)   {
    
    k = kIn;
    scoreList = scoreListIn;
  
    initValue = initValueIn;
    
    bedFileName = bedFileIn;
    min = minIn;
	max = maxIn;
	ltrMin = ltrMinIn;
	
    minPlateauLen = plateauLenIn;
	diffThresh = gapTolIn;
	gapTol = gapTolIn;
	identity = id;

    bList = new vector<BackwardTr *>();
    
    cleanAndMerge();


} 

void MatchTr::cleanAndMerge(){
    
	int len = scoreList->size();
    //cout<<"len="<<len<<endl;
	//cout<<"INITIAL SCORE"<<initValue<<endl;
  
	for(int i = 0;i<len;i++){
        
	    int curr = scoreList->at(i);
        int next;
		
		if(curr != initValue){
			
			int peakLength = 1;
			
			for(int j = i+1;j<len;j++){  
				next = scoreList->at(j);
				if(next == initValue){
					break;
				}
			    peakLength++;
				
			}
			if(peakLength < minPlateauLen){
				
				std::tuple<char,int, int> discard = make_tuple('D',i, i + peakLength);
				//cout<<"Delete pushback"<<endl;
				//cout<<"("<<std::get<0>(discard)<<","<<std::get<1>(discard)<<")"<<endl;
				spikes.push_back(discard);

			}
			else{
				std::tuple<char,int, int> keep = make_tuple('K',i, i + peakLength);
				//cout<<"Keep pushback"<<endl;
               // cout<<"("<<std::get<0>(keep)<<","<<std::get<1>(keep)<<","<<std::get<2>(keep)+k<<")"<<endl;
                spikes.push_back(keep);
			}
			
			i += peakLength;
		}
	}
	//cout<<"Passed initial collections!"<<endl;
	//cout<<"Spikes size:"<<spikes.size()<<endl;
	forwardMerge();
	//cout<<"passed forward Merge"<<endl;

	backwardMerge();
	//cout << "passed bmatchesackward Merge" << endl;
	//medianSmooth();
	smoothAndCollect();
	//cout<<"passed smooth and collect"<<endl;

}


void MatchTr::forwardMerge(){
    //int GAP_TOLERANCE = k;int init = scorer->getInitialScore();
    
    for (int i = 0; i < spikes.size()-1; i++)
	{
		char curr_type;
		int curr_start;
		int curr_end;
		std::tie(curr_type, curr_start, curr_end) = spikes.at(i);


			//	vector<int> section(scoreList->begin()+curr_start,scoreList->begin()+curr_end);

		int level = findMedian(curr_start, curr_end);

		//cout <<"CURR== "<< curr_type << "," << curr_start << "," << curr_end << "," << level << endl;

		char next_type;
        int next_start;
		int next_end;
		std::tie(next_type, next_start, next_end) = spikes.at(i + 1);

		
		//	vector<int> section1(scoreList->begin() + next_start, scoreList->begin() + next_end);
		
	
		int neighborScore = findMedian(next_start,next_end);

	   // cout<<"NEXT= "<<next_type<<","<<next_start<<","<<next_end<<","<<neighborScore<<endl;

		//cout << "level:" << level << ", neighbor:" << neighborScore << endl;
		if (next_type == 'K')
		{
			if (abs(neighborScore - level) < diffThresh && (curr_end + gapTol) >= next_start) // spikes are at same level and within distance of gapTol
			{   
				
				for (int j = curr_end; j <= next_start; j++) //replace curr_end with curr_start
				{
					
					(*scoreList)[j] = neighborScore;
				}
				if (curr_type == 'D')
				{
					spikes[i] = make_tuple('K', curr_start,next_start-1); //flip curr to a section to keep

			
				}
			}
/*
			else if(curr_type == 'K' && (curr_end+gapTol)>=next_start)// two extended sections not at the same level but withing distance of gapTol
			{
				for(int j = curr_end;j<=next_end;j++){  //overwrite old keep section
					(*scoreList)[j] = neighborScore;
				}
			}

*/
				
		
		}

			else if (curr_type =='K' && next_type == 'D'){
			 
				if (abs(neighborScore - level) < diffThresh && (curr_end + gapTol) >= next_start) // spikes are at same level and within distance of k

				{

					
				//cout << "K->D" << endl;

					//  cout<<"MERGING!"<<endl;
					for (int j = curr_end; j <= next_start; j++) //replaced curr_end with curr_start. TODO: Find out why next-start bounds are off at 1824000-1832000
					{
						
						(*scoreList)[j] = level;
					}
					
						spikes[i + 1] = make_tuple('K', curr_start, next_end); //flip next to a section to keep
				}
				
			}
			
	}

}

int MatchTr::findMedian(int start, int end){
	vector<int> section1(scoreList->begin() + start, scoreList->begin() + end);
	if(section1.size() ==0){
		return 0;
	}
	else if(section1.size() ==1){
		return section1.at(0);
	}
	else if (section1.size()%2 ==0){
	std::nth_element(section1.begin(), section1.begin() + section1.size() / 2, section1.end());
	
	return section1[section1.size() / 2];
	}
	else{
		std::nth_element(section1.begin(), section1.begin() + section1.size() / 2, section1.end());

		return section1[(section1.size() / 2)+1];
	}

}

void MatchTr::medianSmooth()
{

	vector<int> *temp = new vector<int>();

	int step = 20;
	int i = 0;

	while (i < scoreList->size() - step)
	{

		int median = findMedian(i, i + step);

		for (int j = i; j < i + step; j++)
		{
			temp->push_back(median);
		}

		i += step;
	}
	int score = findMedian(i, scoreList->size());
	for (int j = i; j < scoreList->size(); j++)
	{
		temp->push_back(score);
	}

	scoreList->clear();
	scoreList = temp;
}

void MatchTr::backwardMerge()
{

	for (int i = spikes.size()-1; i>=1; i--)
	{
		char curr_type;
		int curr_start;
		int curr_end;
		std::tie(curr_type, curr_start, curr_end) = spikes.at(i);

		//vector<int> section(scoreList->begin() + curr_start, scoreList->begin() + curr_end);
	
		//std::nth_element(section.begin(), section.begin() + section.size() / 2, section.end());

		int level = findMedian(curr_start, curr_end);
		

		char next_type;
        int next_start;
		int next_end;
		std::tie(next_type, next_start, next_end) = spikes.at(i - 1);

		//vector<int> section1(scoreList->begin() + next_start, scoreList->begin() + next_end);
		//std::nth_element(section1.begin(), section1.begin() + section1.size() / 2, section1.end());
		int neighborScore = findMedian(next_start,next_end);
		//cout<<"level:"<<level<<", neighbor:"<<neighborScore<<endl;

		if (curr_type == 'K' && next_type == 'D'){

			if (abs(neighborScore - level) < diffThresh && (curr_start - gapTol) <= next_end) // spikes are at same level and within distance of k
			{
				//  cout<<"MERGING!"<<endl;
				for (int j = curr_start; j > next_start; j--)
				{
					(*scoreList)[j] = level;
				}
			}

		}
	}

	   //delete spikes;
}

void MatchTr::smoothAndCollect(){
	//Repeat collection
	int len = scoreList->size();
	//int len = 16000000;
	//int lencandidate = 1000000;
	//vector<std::tuple<int, int>> plateaus;

	for (int i = 0; i < len; i++)
	{

		int curr = scoreList->at(i);
		int next;

		if (curr != initValue)
		{

			int peakLength = 1; 

			for (int j = i + 1; j < len; j++)
			{
				next = scoreList->at(j);
			
				if (next == initValue)
				{
					break;
				}
				peakLength++;
			}
		
			int minSizeKeep = identity * ltrMin / 100;


			//cout << "Peak length is:" << peakLength << endl;
			if (peakLength >= minSizeKeep) //added this parameter
			{
                int height = findMedian(i,i+peakLength-1);

				

				Candidate *keep = new Candidate(i, i + peakLength - 1, height);
				//cout << "(" << std::get<0>(keep) << "," << std::get<1>(keep) << ")" << endl;
				plateaus.push_back(keep);

			}
			else
			{

				for (int k = i; k < i + peakLength; k++)
				{
					(*scoreList)[k] = initValue;
				}

			}
			i += peakLength - 1;
		}
	}
	//	cout<<"SIZE OF PLATEAUS VECTOR"<<plateaus.size()<<endl;
		//std::stack<std::tuple<int,int>> matches; 
       
       PairContainer * matcher = new PairContainer(min,max,diffThresh);
		
		
       // vector<ForwardTr *> * nList = new vector<ForwardTr *>*;		
        Candidate * curr;
		Candidate * match;
		BackwardTr * pair;
		for (int i = 0; i < plateaus.size(); i++)
		{
			curr = plateaus.at(i);
		
			//cout<<"curr:"<<curr->getStart()<<","<<curr->getEnd()<<","<<curr->getHeight()<<endl;
			
			//int currStart = curr.getStart();
			//int currEnd = curr.getEnd();
			//int currHeight = curr.getHeight();
           // cout<<"currStart:"<<currStart<<"currEnd"<<currEnd<<endl;

			/*if(matches.size()>0){
				
				int stackStart;
				int stackEnd;
				std::tie(stackStart, stackEnd) = matches.top();
				//cout<<stackStart<<","<<stackEnd<<endl;

				if (isMatch(stackStart,currStart))
				{
					matches.pop();
					
					ForwardTr *match = new ForwardTr(stackStart, stackEnd, currStart, currEnd);
					fList->push_back(match);
                    
					
				}
				elsecandidate
				{
					matches.push(curr);
					
				}
			}
			else{
				matches.push(curr);
				
			}*/
		//	cout<<"Entering hash"<<endl;
             match = matcher->hashOrReturn(curr);
		//	cout<<"Returning from hash"<<endl;

			if(match!=nullptr){
			//	cout<<"Inside match found"<<endl;
				//cout<<"s1->"<<match->getStart()<<"e1->"<<match->getEnd()<<"s2->"<<curr->getStart()<<"e2->"<<curr->getEnd()<<endl;
				pair = new BackwardTr(match->getStart(), match->getEnd(), curr->getStart(), curr->getEnd());
				//cout<<pair->toString()<<endl;
				bList->push_back(pair);
				//delete pair;
			
			}
			//delete curr;
			//delete match;
	}
	//matcher->empty();

		
} 
bool MatchTr::isMatch(int firstStart,int secondStart){

	int firstHeight = scoreList->at(firstStart);


	int matchLoc = firstStart + firstHeight;
  //  cout<<"matchLoc:"<<matchLoc<<" secondStart:"<<secondStart<<endl;
	return matchLoc == secondStart;

}

vector<BackwardTr *> * MatchTr::getRepeatCandidates(){
    return bList;
}

MatchTr::~MatchTr(){
scoreList->clear();
//delete scoreList;

bList->clear();
delete bList;

}

void MatchTr::printFinalScores(int start, int end)
{

	ofstream output;

	output.open("../output/cleanedScores.csv");
	int len = end > scoreList->size() ? scoreList->size() : end;

	for (int i = start; i < len; i++)
	{

		output << i << "," << scoreList->at(i) << endl;
	}
	output.close();
}

void MatchTr::bedFormat(int start, int end)
{
    cout<<bedFileName<<endl;
    ofstream output;
    output.open(bedFileName);

    for (int i = 0; i < bList->size(); i++)
    {
        BackwardTr *curr = bList->at(i);

        output << "chrX"<< "	" << curr->getS1() << "	" << curr->getE2()<< endl;
    }
    output.close();
}
}