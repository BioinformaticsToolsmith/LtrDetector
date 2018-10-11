#include <iostream>
#include <TrCollector.h>

using namespace std;
using namespace nonltr;
using namespace utility;
using namespace exception;

int main(){

    ChromosomeOneDigit * chrX = new ChromosomeOneDigit("~/Downloads/Dm6/Test");

    int kval = 10;

    int min = ?;
    int max =?;
    int d = ?;

    


    TrCollector * app = new TrCollector(chrX, kval, min, max, d);

    

    return 0;
}