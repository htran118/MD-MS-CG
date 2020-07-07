#include <iostream>
#include <vector>
#include <stdlib.h>
#include "../libs/const.h"

using namespace MD;
using namespace std;

int main(int argc, char *argv[]) {
  vector<myVect*> rCor;
  vector<myCoor> fCor;
  int mSze, leng, pNum;
  myCoor rand, fric;
  vector<int> bins, dBin;

  if(argc != 4) {
    cerr << "Usage:" << argv[0] << " <leng> <mSze> <pNum>" << endl;
    return 1;
  }

  leng = atoi(argv[1]);
  mSze = atoi(argv[2]);
  pNum = atoi(argv[3]);
  pNum = pNum * (pNum - 1) / 2;

  for(int i = 0; i < (leng * pNum); ++i) {
    rCor.push_back(myVectMalloc(DIM));
    for(int k = 0; k < DIM; ++k) {
      cin >> rand;
      myVectSet(rCor.at(i), k, rand);
    }
    //myVectSet(rCor.at(i), 0, rand);
  }

  for(int i = 0; i < mSze; ++i) {
    fric = 0.0;
    for(int k = 0; k < pNum; ++k)
      for(int j = 0; j < (leng - i); ++j) {
        myVectDot(rCor.at(j * pNum + k), 
                  rCor.at((i + j) * pNum + k), &rand);
        fric += rand;
        //fric += myVectGet(rCor.at(j * pNum + k), 0) * 
        //        myVectGet(rCor.at((i + j) * pNum + k), 0);
    }
    cout << (fric / (leng - i) / pNum * STPSZE) << endl;
  }
  
  return 0;
}
