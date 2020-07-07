#include <iostream>
#include <stdlib.h>
#include "../libs/const.h"

using namespace std;
using namespace MD;

int main() {
  myVect *vec1 = myVectMalloc(DIM), *vec2 = myVectMalloc(DIM);
  myCoor sdot;
  normal_distribution<myCoor> norm{0.0,1.0};
  char again = 'y';
  while(again == 'y') {
    for(int i = 0; i < DIM; ++i) {
      myVectSet(vec1, i, norm(ENGINE));
      myVectSet(vec2, i, norm(ENGINE));
    }
    myVectOut(vec1);
    myVectOut(vec2);
    myVectCrs(vec1, vec2);
    myVectOut(vec1);
    myVectDot(vec1, vec2, &sdot);
    cout << sdot << endl;
    cout << "Again? (y/n)" << endl;
    cin >> again;
  }
  return 0;
}
