#include <cmath>
#include <cassert>
#include <stdlib.h>
#include "const.h"

using namespace std;

namespace MD
{
  extern const myCoor ONEPIE = acos((double) -1.0);
  extern const myCoor TWOPIE = ONEPIE * 2;
  extern const myCoor ROOTWO = 1.41421;
  extern const myCoor IVRTWO = 1 / ROOTWO;
  extern const myCoor COULMB = 1;
  extern const myCoor BOLTZM = 0.00831451;
  extern const myCoor MINZRO = pow(10,-5);

  extern const myCoor STPSZE = 0.001000;
  extern const myCoor STPHLF = STPSZE * 0.5;
  extern const myCoor STPINV = 1 / STPSZE;
  extern const myCoor STPIHF = 1 / STPHLF;
  extern const myCoor STEPSQ = STPSZE * STPSZE;
  extern const myCoor STPIRT = 1 / sqrt(STPSZE);
  extern const int DIM = 3;

  extern const int MAXNPART = 120;
  extern const int MAXNBOND = MAXNPART * 2;
  extern const int MAXNANGL = MAXNBOND * 2;
  extern const int MAXNDHED = MAXNBOND * 2;

  extern const int BONDFORC = 2; 
  extern const bool ANGLPART = false;
  extern const bool LJONPART = true;
  extern const bool DPDPARTK = true;

  extern const char *defNme = "define/";
  extern const char *datNme = "data/";
  extern const char *paraFl = "define/param.txt";
  extern const char *topoFl = "define/topol.txt";
  extern const char *strcFl = "define/struc.txt";
  extern const char *cgStFl = "define/cgStr.txt";
  extern const char *coorFl = "define/coord.txt";
  extern const char *trajFl = "data/trajc.trj";
  extern const char *memoFl = "define/memor.txt";
  extern const char *convFl = "data/coord.crd";

  extern const char *crdFmt = "%9.4f";
  extern const char *typFmt = "%-5c";
  extern const char *sNmFmt = "%-5s";
  extern const char *nmeFmt = "%-5d";

  extern const string ctrlAtom = "ATOM";
  extern const string ctrlMass = "MASS";
  extern const string ctrlBond = "BOND";
  extern const string ctrlAngl = "ANGL";
  extern const string ctrlDhed = "DHED";
  extern const string ctrlRsid = "RSID";
  extern const string ctrlTime = "TIME";
  extern const string ctrlStop = "STOP";
  extern const string ctrlIgnr = "!";
  extern const string ctrlEnd  = "END";
  extern const string emptyNme = "NaN";

  void myVectInv(myVect* vec1, const myCoor leng)
  { myVectScl(vec1, 1 / leng); }

  void myVectRev(myVect* vec1)
  { myVectScl(vec1, -1); }

  //get the cross product of 2 vectors, with exception handle
  //for parallel case
  void myVectCrs(myVect* vec1, const myVect* vec2) {
    assert(DIM == 3);
    myCoor x = myVectGet(vec1, 1) * myVectGet(vec2, 2) -
               myVectGet(vec1, 2) * myVectGet(vec2, 1),
           y = myVectGet(vec1, 2) * myVectGet(vec2, 0) -
               myVectGet(vec1, 0) * myVectGet(vec2, 2),
           z = myVectGet(vec1, 0) * myVectGet(vec2, 1) -
               myVectGet(vec1, 1) * myVectGet(vec2, 0);
    if(abs(x) < MINZRO && abs(y) < MINZRO && abs(z) < MINZRO)
      myVectOrd(vec1, vec2);
    else {
      myVectSet(vec1, 0, x);
      myVectSet(vec1, 1, y);
      myVectSet(vec1, 2, z);
    }
  }

  void myVectVecCrs(const myVect* vec1, const myVect* vec2, myVect* vec3) {
    assert(DIM == 3 && vec1 != vec3 && vec2 != vec3);
    myVectSet(vec3, 0, (myVectGet(vec1, 1) * myVectGet(vec2, 2) -
                        myVectGet(vec1, 2) * myVectGet(vec2, 1)));
    myVectSet(vec3, 1, (myVectGet(vec1, 2) * myVectGet(vec2, 0) -
                        myVectGet(vec1, 0) * myVectGet(vec2, 2)));
    myVectSet(vec3, 2, (myVectGet(vec1, 0) * myVectGet(vec2, 1) -
                        myVectGet(vec1, 1) * myVectGet(vec2, 0)));
    if(abs(myVectGet(vec3, 0)) < MINZRO &&
       abs(myVectGet(vec3, 1)) < MINZRO &&
       abs(myVectGet(vec3, 2)) < MINZRO)
      myVectOrd(vec3, vec1);
  }

  //scale to get unit vector, with exception handle for zero value
  //POSTCON: dist is the original length of the vector
  void myVectUnt(myVect* uVec, myCoor& dist) {
    myCoor leng = myVectNrm(uVec);
    dist = leng;
    while(leng < MINZRO) {
      myVectRnd(uVec);
      leng = myVectNrm(uVec);
    }
    myVectInv(uVec, leng);
  }

  //get a random vector orthogonal to another
  void myVectOrd(myVect* uVec, const myVect* sVec) {
    myCoor dist;
    switch(DIM) {
      case 1:
        myVectZro(uVec);
        break;
      case 2:
        dist = NRMDIS(ENGINE);
        myVectSet(uVec, 0, dist * myVectGet(sVec, 1));
        myVectSet(uVec, 1, -dist * myVectGet(sVec, 0));
        break;
      case 3:
        for(int k = 0; k < DIM; ++k)
          myVectSet(uVec, k, NRMDIS(ENGINE));
        myVectCrs(uVec, sVec);
        break;
      default:
        cerr << "Not supporting current dimension " << DIM << endl;
        break;
    }
  }

  //generate a random vector
  void myVectRnd(myVect* uVec) {
    for(int i = 0; i < (*uVec).size; ++i)
      myVectSet(uVec, i, NRMDIS(ENGINE)); 
  }

  void myVectInp(myVect* vect, istream& iput) {
    myCoor dumy;
    for(int i = 0; i < (*vect).size; ++i) {
      iput >> dumy;
      myVectSet(vect, i, dumy);
    }
  }

  void myVectOut(const myVect* vect) {
    for(int i = 0; i < (*vect).size; ++i)
      COUTPT << " " << myVectGet(vect, i);
    COUTPT << " ";
  }

  void dbuggr() {
    cout << "FUCK" << endl;
  }
}
