#include <iostream>
#include <vector>
#include <stdlib.h>
#include "../libs/const.h"

using namespace MD;
using namespace std;

int main(int argc, char *argv[]) {
  vector <myVect*> rand;
  myCoor vdot, iput;
  int memL, memO;

  if(argc != 3) {
    cerr << "Usage:" << argv[0] << " <memLen> <memOff>" << endl;
    return 1;
  }

  memL = atoi(argv[1]);
  memO = atoi(argv[2]);

  for(int i = 0; i < memL; ++i) {
    rand.push_back(myVectMalloc(DIM));
    for(int k = 0; k < DIM; ++k) {
      cin >> iput;
      myVectSet(rand.at(i), k, iput);
    }
  }

  for(int i = 0; i < memO; ++i) {
    vdot = 0.0;
    for(int j = 0; j < (memL - i); ++j) {
      myVectDot(rand.at(j), rand.at(i + j), &iput);
      vdot += iput;
    }
    vdot *= STPSZE / (memL - i);
    cout << vdot << endl;
  }

  return 0;
}
