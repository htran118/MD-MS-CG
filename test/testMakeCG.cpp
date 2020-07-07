/****************************************************************************
 * Convert trajectory files to standard ouput
 * Default input file is data/trajc.trj
 ***************************************************************************/

#include <iostream>
#include <sstream>
#include <string>
#include <cstring>
#include <fstream>
#include <cmath>
#include "../libs/const.h"

using namespace std;
using namespace MD;

const int RSIDNO = 125;
const int RSIDSZ = 40;
const int CGSIZE = 10;
const myCoor MB = 15.035;
const myCoor MA = 12.027;

int main(int argc, char *argv[]) {
  char optFle[256], outFle[256], inpFle[256], format[256];
  myCoor dummy, mass, mTot;
  myVect *xTot = myVectMalloc(DIM), *vTot = myVectMalloc(DIM),
         *dumy = myVectMalloc(DIM);
  string dump;
  FILE *pf;

  if(argc != 3) {
    cerr << "Usage:" << argv[0] << " <inpt> <outp>" << endl;
    cerr << "Using default path" << endl;
    strcpy(inpFle, coorFl);
    strcat(inpFle, ".gro");
    strcpy(outFle, coorFl);
  }
  else {
    strcpy(inpFle, argv[1]);
    strcpy(outFle, argv[2]);
  }
  cout << argv[0] << " < " << inpFle << " > " << outFle << endl;
    
  ifstream inpt(inpFle, ios_base::in);
  pf = fopen(outFle, "w");

  inpt.ignore(256, '\n');
  inpt.ignore(256, '\n');
  fprintf(pf, "0\n");
  
  for(int i = 0; i < RSIDSZ * RSIDNO / CGSIZE; ++i) {
    mTot = 0;
    myVectZro(xTot);
    myVectZro(vTot);
    for(int j = 0; j < CGSIZE; ++j) {
      if((i * CGSIZE + j) % RSIDSZ == RSIDSZ - 1||
         (i * CGSIZE + j) % RSIDSZ == 0)
        mass = MB;
      else
        mass = MA;
      mTot += mass;
      for(int k = 0; k < 3; ++k)
        inpt >> dump;
      for(int k = 0; k < DIM; ++k) {
        inpt >> dummy;
        myVectSet(dumy, k, dummy);
      }
      myVectScl(dumy, mass);
      myVectAdd(xTot, dumy);
      for(int k = 0; k < DIM; ++k) {
        inpt >> dummy;
        myVectSet(dumy, k, dummy);
      }
      myVectScl(dumy, mass);
      myVectAdd(vTot, dumy);
    }
    myVectInv(xTot, mTot);
    myVectInv(vTot, mTot);
    for(int k = 0; k < DIM; ++k)
      fprintf(pf, crdFmt, myVectGet(xTot, k));
    for(int k = 0; k < DIM; ++k)
      fprintf(pf, crdFmt, myVectGet(vTot, k));
    fprintf(pf, "\n");
  }
  return 0;
}
