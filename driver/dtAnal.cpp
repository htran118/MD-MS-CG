/****************************************************************************
 * Data analyzer
 * Default input file is data/trajc.trj
 ***************************************************************************/

#include <iostream>
#include <sstream>
#include <string>
#include <cstring>
#include <fstream>
#include <cmath>
#include "../libs/const.h"
#include "../libs/tlist.h"
#include "../libs/world.h"

using namespace MD;
using namespace std;

int PBINSZ = 500;
int ABINSZ = 200;
int BBINSZ = 200;

int main(int argc, char *argv[]) {
  char optTrj[256], optOut[256];
  world systm;
  tlist types;
  vector<partk> *parts;
  vector<myVect*> oPos, oVel, oPVl; 
  int duration, interval, sum;
  myCoor vacf, msdf, pvac, kErg, pErg, ctOf;
  vector<int> pDum, bDum, aDum, dDum;
  vector<vector<int>> pBin, bBin, aBin, dBin;
  myVect* vDum;

  if(argc != 2) {
    cerr << "Usage:" << argv[0] << " <pathname>" << endl;
    cerr << "Using default path" << endl;
    strcpy(optTrj, trajFl);
  }
  else {
    strcpy(optTrj, argv[1]);
  }
 
  COUTPT << "Reading parameter file..." << endl;
  ifstream parameter(paraFl, ios_base::in);
  types.adType(parameter);
  parameter.close();
  COUTPT << "Reading topology file..." << endl;
  ifstream topology(topoFl, ios_base::in);
  types.adType(topology);
  topology.close();
  COUTPT << "Reading structure file..." << endl;
  ifstream structure(strcFl, ios_base::in);
  systm.iniSys(structure, &types);
  structure.close();

  COUTPT << "Initialization..." << endl;  
  parts = (const_cast<vector<partk>*> (&(systm.obPartLs())));
  for(int i = 0; i < systm.obNpartk(); ++i) {
    oPos.push_back(myVectMalloc(DIM));
    oVel.push_back(myVectMalloc(DIM));
    for(int j = i + 1; j < systm.obNpartk(); ++j) {
      oPVl.push_back(myVectMalloc(DIM));
    }
  }
  for(int i = 0; i < PBINSZ; ++i)
    pDum.push_back(0);
  for(int i = 0; i < types.obNptype(); ++i)
    for(int j = 0; j < types.obNptype(); ++j) {
      pBin.push_back(vector<int> (0));
      for(int l = 0; l < PBINSZ; ++l)
        (pBin.at(i * types.obNptype() + j)).push_back(0);
    }
  for(int i = 0; i < BBINSZ; ++i)
    bDum.push_back(0);
  for(int i = 0; i < types.obNbtype(); ++i) {
    bBin.push_back(vector<int> (0));
    for(int j = 0; j < BBINSZ; ++j)
      (bBin.at(i)).push_back(0);
  }
  for(int i = 0; i < ABINSZ; ++i)
    aDum.push_back(0);
  for(int i = 0; i < types.obNatype(); ++i) {
    aBin.push_back(vector<int> (0));
    for(int j = 0; j < ABINSZ; ++j)
      (aBin.at(i)).push_back(0);
  }

  COUTPT << "Analysing trajectory..." << endl;
  ifstream trajectory(optTrj, ios_base::in);
  systm.iniCon(trajectory);
  systm.mkDcon(0);
  trajectory >> duration >> interval;
  cout << duration << " " << interval << endl;
  for(int n = 0; n < (duration / interval); ++n) {
    systm.iniTrj(trajectory);
    cout << "Reading frame at t=" << systm.obTime() << endl;
    if(n == 0)
      for(int i = 0; i < systm.obNpartk(); ++i) {
        myVectCpy(oPos.at(i), ((*parts).at(i)).obPost());
        myVectCpy(oVel.at(i), ((*parts).at(i)).obVelo());
        //for(int j = i + 1; j < systm.obNpartk(); ++j) {
        //  obVdif(&((*parts).at(j)), &((*parts).at(i)),
        //         oPVl.at(obPInd(parts, &((*parts).at(j)), &((*parts).at(i)))));
        //}
      }
    msdf = systm.obtMSD(&oPos);
    vacf = systm.obtVAC(&oVel);
    //pvac = systm.obPVVC(&oPVl);
    kErg = systm.obKine();
    //pErg = systm.obPote();
    for(int i = 0; i < types.obNptype(); ++i)
      for(int j = 0; j < types.obNptype(); ++j) {
        pDum = systm.obLjDs(types.obPtyp(i), types.obPtyp(j), PBINSZ);
        for(int m = 0; m < PBINSZ; ++m)
          (pBin.at(i * types.obNptype() + j)).at(m) += pDum.at(m);
      }
    for(int i = 0; i < types.obNbtype(); ++i) {
      bDum = systm.ob12Ds(types.obBtyp(i), BBINSZ);
      for(int m = 0; m < BBINSZ; ++m)
        (bBin.at(i)).at(m) += bDum.at(m);
    }
    for(int i = 0; i < types.obNatype(); ++i) {
      aDum = systm.ob13Ds(types.obAtyp(i), ABINSZ);
      for(int m = 0; m < ABINSZ; ++m)
        (aBin.at(i)).at(m) += aDum.at(m);
    }
    cout << n * interval * STPSZE 
         << '\t' << msdf
         << '\t' << vacf
         //<< '\t' << pvac
         << '\t' << kErg
         //<< '\t' << pErg
         //<< '\t' << kErg - pErg
         << '\t' << (*(systm.obBond(0))).obSize() 
         << endl;
  }
  ctOf = max(systm.obOfLj(), systm.obOfDp());
  ctOf = (ctOf < MINZRO) ? (myVectMin(systm.obEdge()) / 2) : (ctOf * 2);
  for(int i = 0; i < types.obNptype(); ++i)
    for(int j = 0; j < types.obNptype(); ++j) {
      sum = 0;
      for(int m = 0; m < PBINSZ; ++m)
        sum += (pBin.at(i * types.obNptype() + j)).at(m);
      if(sum > 0) {
        cout << endl << "Distribution for non-bonded interaction between "
             << (*(types.obPtyp(i))).obName() << " and "
             << (*(types.obPtyp(j))).obName() << endl;
        for(int m = 0; m < PBINSZ; ++m)
          cout << ctOf * (m + 0.5) / PBINSZ << '\t'
               << (pBin.at(i * types.obNptype() + j)).at(m) << endl;
      }
      else
        if(systm.obLcon() != 0)
          cout << endl << "No non-bonded interaction between "
               << (*(types.obPtyp(i))).obName() << " and "
               << (*(types.obPtyp(j))).obName() << endl;
    }
  for(int i = 0; i < types.obNbtype(); ++i) {
    const btype *btyp = types.obBtyp(i);
    sum = 0;
    for(int j = 0; j < BBINSZ; ++j)
      sum += (bBin.at(i)).at(j);
    if(sum > 0) {
      cout << endl << "Distribution for 1-2 interaction between "
           << (*((*btyp).obLeft())).obName() << " and "
           << (*((*btyp).obRigh())).obName() << endl;
      for(int j = 0; j < BBINSZ; ++j)
        cout << 2 * (*btyp).obRelx() * (j + 0.5) / BBINSZ
             << '\t' << (bBin.at(i)).at(j) << endl;
    }
    else
      if(BONDFORC > 0)
        cout << endl << "No 1-2 interaction between "
             << (*((*btyp).obLeft())).obName() << " and "
             << (*((*btyp).obRigh())).obName() << endl;
  }
  for(int i = 0; i < types.obNatype(); ++i) {
    const atype *atyp = types.obAtyp(i);
    sum = 0;
    for(int j = 0; j < ABINSZ; ++j)
      sum += (aBin.at(i)).at(j);
    if(sum > 0) {
      cout << endl << "Distribution for 1-3 interaction between "
           << (*((*atyp).obMidl())).obName() << " to "
           << (*((*atyp).obRigh())).obName() << " and "
           << (*((*atyp).obLeft())).obName() << endl;
      for(int j = 0; j < ABINSZ; ++j)
        cout << (1.0 * j / ABINSZ - 0.5) * 2
             << '\t' << (aBin.at(i)).at(j) << endl;
    }
    else
      if(BONDFORC > 1)
        cout << endl << "No 1-3 interaction between "
             << (*((*atyp).obMidl())).obName() << " to "
             << (*((*atyp).obRigh())).obName() << " and "
             << (*((*atyp).obLeft())).obName() << endl;
  }
  return 0;
}
