#include <iostream>

#include "PairContainer.h"
#include "Candidate.h"
using namespace std;
using namespace tr;

int main(){
    cout<<"Building matcher"<<endl;
    PairContainer * matcher = new PairContainer(2000,36000,150);
    cout<<"Matcher made"<<endl;
    Candidate  * first = new Candidate(8000,9000,13850);
    cout<<"First Candidate made"<<endl;
    Candidate  * second = new Candidate(57195,57609,6768);
    cout<<"Second Candidate made"<<endl;
    Candidate * third = new Candidate(63963,64377,-6768);

    Candidate * fourth = new Candidate(17000,19000,13800);

    Candidate * firstAttempt = matcher->hashOrReturn(first);
    cout<<"First hash passed"<<endl;
    if(firstAttempt ==nullptr){
        cout<<"Great Success!"<<endl;
    }
    Candidate * secondAttempt = matcher->hashOrReturn(second);
    cout<<"Second hash passed"<<endl;
   // cout<<secondAttempt->getStart()<<secondAttempt->getEnd()<<secondAttempt->getHeight()<<endl;

    Candidate * thirdAttempt = matcher->hashOrReturn(third);

        cout << thirdAttempt->getStart() << thirdAttempt->getEnd() << thirdAttempt->getHeight() << endl;
    

    Candidate * fourthAttempt = matcher->hashOrReturn(fourth);
        cout<<fourthAttempt->getStart()<<fourthAttempt->getEnd()<<fourthAttempt->getHeight()<<endl;
}
