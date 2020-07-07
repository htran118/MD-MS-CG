/********************************************************************
 * Memory kernel extractor
 *******************************************************************/

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

myCoor RCTOFF = 1.00;
myCoor RSTPSZ = 0.02;
int MSKPSZ = 50;
int MCTOFF = 50; 
int MRSIZE = (int) (RCTOFF / RSTPSZ);
myVect *FORM = myVectMalloc(DIM), *UNIM = myVectMalloc(DIM),
       *TOTM = myVectMalloc(DIM);

void obDelF(world *aSyst, world *cgSys, const int i, const int j,
            const int cgRt, myVect* delF, myCoor& para, myCoor& orth);

int main(int argc, char *argv[]) {
  char optTrj[256], optOut[256];
  world aSyst, cgSys;
  tlist types;
  partk *atom1, *atom2, *part1, *part2; 
  vector<resid> *rsids;
  vector<partk> *aPart, *cgPrk;
  vector<myCoor> oDis;
  vector<bool> pOut;
  vector<myVect*> oPos, oVel, delF;
  vector<myCoor> pCorNbnd, pCorBond, oCorNbnd, oCorBond, oFor, pFor;
  vector<int> contNbnd, contBond;
  myVect *vDum, *xDum, *fDum;
  myMatr *mDum;
  myCoor nDis, para, orth, dist, pDum, oDum;
  int duration, interval, sum, trjC, sysC, cgRt, cgRsidNo, partPRsd, rLoc,
      indx;
  FILE *pf;

  if(argc != 2) {
    cerr << "Usage:" << argv[0] << " <pathname>" << endl;
    cerr << "Using default path" << endl;
    strcpy(optTrj, trajFl);
  }
  else {
    strcpy(optTrj, argv[1]);
  }

  COUTPT << "This program extracts the memory kernel from an atomistic" << endl
         << "trajectory. Currently support POL40" << endl
         << endl;
  COUTPT << "The memory cutoff is " << MCTOFF << " and the radius cutoff is "
         << MRSIZE << endl;

  COUTPT << "Choose the coarse-graining rate" << endl;
  CINPUT >> cgRt;

  COUTPT << "Reading parameter file..." << endl;
  ifstream parameter(paraFl, ios_base::in);
  types.adType(parameter);
  parameter.close();
  COUTPT << "Reading topology file..." << endl;
  ifstream topology(topoFl, ios_base::in);
  types.adType(topology);
  topology.close();
  if((*(types.obRtyp("POL40"))).obNptype() % cgRt != 0) {
    cerr << "Coarse-graining rate does not divide length of POL40" << endl;
    terminate();
  }
  COUTPT << "Reading structure file..." << endl;
  ifstream structure(strcFl, ios_base::in);
  aSyst.iniSys(structure, &types);
  structure.close();

  partPRsd = (*(types.obRtyp("POL40"))).obNptype() / cgRt;
  cgRsidNo = aSyst.obNresid();
  cout << cgRsidNo << " " << partPRsd << endl;
  COUTPT << "Generating CG structure file..." << endl;
  pf = fopen(cgStFl, "w");
  for(int i = 0; i < cgRsidNo; ++i)
    fprintf(pf, "%s %-5d%s\n", ctrlRsid.c_str(), (i + 1), "CPN");
  fprintf(pf, "%s\n", ctrlStop.c_str());
  fprintf(pf, "\n%s List of particles\n", ctrlIgnr.c_str());
  fprintf(pf, "%s    NAME RNAM TYPE\n", ctrlIgnr.c_str());
  for(int i = 0; i < cgRsidNo; ++i)
    for(int j = 0; j < partPRsd; ++j)
      fprintf(pf, "%s %-5d%-5d%s%d   \n", ctrlAtom.c_str(),
                                          (partPRsd * i + j + 1),
                                          (i + 1), "C", (j + 1));
  fprintf(pf, "%s\n", ctrlStop.c_str());
  fprintf(pf, "\n%s List of bonds\n", ctrlIgnr.c_str());
  fprintf(pf, "%s    LEFT RIGH\n", ctrlIgnr.c_str());
  for(int i = 0; i < cgRsidNo; ++i)
    for(int j = 0; j < partPRsd - 1; ++j)
      fprintf(pf, "%s %-5d%-5d\n", ctrlBond.c_str(),
                                   (partPRsd * i + j + 1),
                                   (partPRsd * i + j + 2));
  fprintf(pf, "%s\n", ctrlStop.c_str());
  fprintf(pf, "\n%s List of angles\n", ctrlIgnr.c_str());
  fprintf(pf, "%s    LEFT MIDL RIGH\n", ctrlIgnr.c_str());
  for(int i = 0; i < cgRsidNo; ++i)
    for(int j = 0; j < partPRsd - 2; ++j)
      fprintf(pf, "%s %-5d%-5d%-5d\n", ctrlAngl.c_str(),
                                       (partPRsd * i + j + 1),
                                       (partPRsd * i + j + 2),
                                       (partPRsd * i + j + 3));
  fprintf(pf, "\n%s List of dihedral angles\n", ctrlIgnr.c_str());
  fprintf(pf, "%s    UPPR LEFT RIGH LOWR\n", ctrlIgnr.c_str());
  for(int i = 0; i < cgRsidNo; ++i)
    for(int j = 0; j < partPRsd - 3; ++j)
      fprintf(pf, "%s %-5d%-5d%-5d%-5d\n", ctrlDhed.c_str(),
                                           (partPRsd * i + j + 1),
                                           (partPRsd * i + j + 2),
                                           (partPRsd * i + j + 3),
                                           (partPRsd * i + j + 4));
  fprintf(pf, "%s\n", ctrlStop.c_str());
  fprintf(pf, "END\n");
  fclose(pf);
  ifstream cgStruct(cgStFl, ios_base::in);
  cgSys.iniSys(cgStruct, &types);
  cgStruct.close();

  ifstream trajectory(optTrj, ios_base::in);
  aSyst.iniCon(trajectory);
  trajectory >> duration >> interval;
  cout << duration << " " << interval << endl;
  cgSys.mkEdge(aSyst.obEdge());
  trajectory.clear();
  trajectory.seekg(0, ios_base::beg);
  cgSys.iniCon(trajectory);

  COUTPT << "Initialization..." << endl;
  aPart = (const_cast<vector<partk>*> (&(aSyst.obPartLs())));
  cgPrk = (const_cast<vector<partk>*> (&(cgSys.obPartLs())));
  for(int i = 0; i < cgSys.obNpartk(); ++i) 
    for(int j = i + 1; j < cgSys.obNpartk(); ++j) {
      delF.push_back(myVectMalloc(DIM));
      oVel.push_back(myVectMalloc(DIM));
      oPos.push_back(myVectMalloc(DIM));
      oFor.push_back(0.0);
      pFor.push_back(0.0);
      oDis.push_back(0.0);
      pOut.push_back(false);
    }
  for(int m = 0; m < MCTOFF; ++m) 
    for(int n = 0; n < MRSIZE; ++n) {
      pCorNbnd.push_back(0.0);
      pCorBond.push_back(0.0);
      oCorNbnd.push_back(0.0);
      oCorBond.push_back(0.0);
      //corrNbnd.push_back(myMatrMalloc(DIM, DIM));
      //corrBond.push_back(myMatrMalloc(DIM, DIM));
      //myMatrZro(corrNbnd.back());
      //myMatrZro(corrBond.back());
      contNbnd.push_back(0);
      contBond.push_back(0);
    }
  vDum = myVectMalloc(DIM);
  xDum = myVectMalloc(DIM);
  fDum = myVectMalloc(DIM);
  mDum = myMatrMalloc(DIM, DIM);
  trjC = 0;

  COUTPT << "Extracting memory..." << endl;
  for(int n = 0; n < ((duration / interval) - MCTOFF) / ((myCoor) MSKPSZ);
                 ++n) {
    //restart from the beginning
    trajectory.clear();
    trajectory.seekg(0, ios_base::beg);
    aSyst.iniCon(trajectory);
    trajectory >> duration >> interval;
    for(int i = 0; i < cgSys.obNpartk(); ++i) 
      for(int j = i + 1; j < cgSys.obNpartk(); ++j)
        pOut.at(obPInd(cgPrk, &((*cgPrk).at(i)),
                              &((*cgPrk).at(j)))) = false;
    sysC = 0;
    //skip some frames
    for(int m = 0; m < trjC * MSKPSZ; ++m)
      aSyst.skpTrj(trajectory);
    COUTPT << "Origin set at frame " << (trjC * MSKPSZ) << "..." << endl;
    for(int m = 0; m < MCTOFF; ++m) {
      COUTPT << "Reading frame " << m << "..." << endl;
      //if all pairs have been considered, break the loop
      if(sysC >= (cgSys.obNpartk() * (cgSys.obNpartk() - 1) / 2))
        break;
      //input new coordinates for both system
      aSyst.iniTrj(trajectory);
      for(int i = 0; i < (*cgPrk).size(); ++i) {
        myVectZro(xDum);
        myVectZro(vDum);
        for(int j = 0; j < cgRt; ++j) {
          myVectAdd(xDum, ((*aPart).at(cgRt * i + j)).obPost());
          myVectAdd(vDum, ((*aPart).at(cgRt * i + j)).obVelo());
        }
        myVectScl(xDum, 1.0 / cgRt);
        myVectScl(vDum, 1.0 / cgRt);
        ((*cgPrk).at(i)).mkPost(xDum);
        ((*cgPrk).at(i)).mkVelo(vDum);
      }
      //save original coordinates
      if(m == 0)
        for(int i = 0; i < cgSys.obNpartk(); ++i) {
          part1 = &((*cgPrk).at(i));
          for(int j = i + 1; j < cgSys.obNpartk(); ++j) {
            part2 = &((*cgPrk).at(j));
            indx = obPInd(cgPrk, part1, part2);
            cgSys.ljonSz(part1, part2, oDis.at(indx), UNIM);
            if(oDis.at(indx) > RCTOFF) {
              pOut.at(indx) = true;
              ++sysC;
              continue;
            }
            //debugger...
            if(oDis.at(indx) < MINZRO) {
              cerr << "Distance close to zero" << endl;
              cerr << "time=" << aSyst.obTime() << ", pair=" << i << "&" << j
                   << ", dist=" << oDis.at(indx) << endl;
              cerr << "Coordinate=";
              myVectOut((*part1).obPost());
              myVectOut((*part2).obPost());
              cerr << endl;
              terminate();
            }
            obVdif(part1, part2, vDum);
            myVectCpy(oVel.at(indx), vDum);
            obDelF(&aSyst, &cgSys, i, j, cgRt, delF.at(indx), 
                   pFor.at(indx), oFor.at(indx));
          }
        }
      //memory extraction
      for(int i = 0; i < cgSys.obNpartk(); ++i) {
        part1 = &((*cgPrk).at(i));
        for(int j = i + 1; j < cgSys.obNpartk(); ++j) {
          part2 = &((*cgPrk).at(j));
          indx = obPInd(cgPrk, part1, part2);
          if(pOut.at(indx))
            continue;
          //new distance
          cgSys.ljonSz(part1, part2, nDis, UNIM);
          myVectCpy(TOTM, delF.at(indx));
          myVectDot(TOTM, UNIM, &pDum);
          myVectScl(UNIM, pDum);
          myVectSub(TOTM, UNIM);
          oDum = myVectNrm(TOTM);
          //if new distance is too far from old distance, go to next pair
          if(RSTPSZ < abs(nDis - oDis.at(indx))) {
            pOut.at(indx) = true;
            ++sysC;
            continue;
          }
          rLoc = (int) floor(oDis.at(indx) / RSTPSZ);
          //debugger...
          if(rLoc >= MRSIZE) {
            cerr << "rLoc overflown" << endl;
            cerr << "time=" << aSyst.obTime() << ", pair=" << i << "&" << j
                   << ", dist=" << oDis.at(indx) << endl;
            cerr << "Coordinate=";
            myVectOut((*part1).obPost());
            myVectOut((*part2).obPost());
            cerr << endl;
            pOut.at(indx) = true;
            ++sysC;
            continue; 
          }
          //extract memory
          obDelF(&aSyst, &cgSys, i, j, cgRt, fDum, para, orth);
          //debugger...
          if(para > 100000) {
            for(int a = 0; a < cgRt; ++a) {
              atom1 = const_cast<partk*> (aSyst.obPart(i * cgRt + a));
              for(int b = 0; b < cgRt; ++b) {
                pOut.at(indx) = true;
                atom2 = const_cast<partk*> (aSyst.obPart(j * cgRt + b));
                aSyst.ljonSz(atom1, atom2, dist, UNIM);
                aSyst.ljonFc(atom1, atom2, 0, dist, UNIM, FORM); 
                myVectDot(FORM, UNIM, &para);
                if(para > 100000) {
                  cout << i << " " << j << " ";
                  cout << (i * cgRt + a) << " " << (j * cgRt + b) << " ";
                  myVectOut((*atom1).obPost());
                  myVectOut((*atom2).obPost());
                  myVectOut(UNIM);
                  cout << " " << dist << " ";
                  myVectOut(FORM);
                  cout << endl;
                }
              }
            }
          }
          if(!samBnd(part1, part2)) {
              //pCorNbnd.at(m * MRSIZE + rLoc) += para * pFor.at(indx);
              //oCorNbnd.at(m * MRSIZE + rLoc) += orth * oFor.at(indx);
              pCorNbnd.at(m * MRSIZE + rLoc) += para * pDum;
              oCorNbnd.at(m * MRSIZE + rLoc) += orth * oDum;
          }
          else { }
          //for(int k = 0; k < DIM; ++k) {
          //  myVectCpy(xDum, oFor.at(indx));
          //  myVectMul(xDum, fDum);
          //  myMatrSetrow(mDum, k, xDum);
          //}
          ++contNbnd.at(m * MRSIZE + rLoc);
        }
      }
    }
    ++trjC;
  }

  COUTPT << "\n\ncount=" << endl;
  for(int m = 0; m < MCTOFF; ++m) {
    for(int n = 0; n < MRSIZE; ++n)
      COUTPT << contNbnd.at(m * MRSIZE + n) << " ";
   COUTPT << endl; 
  }
  COUTPT << "\n\npCorr=" << endl;
  for(int m = 0; m < MCTOFF; ++m) {
    for(int n = 0; n < MRSIZE; ++n) {
      if(contNbnd.at(m * MRSIZE + n) != 0)
        pCorNbnd.at(m * MRSIZE + n) *= 1.0 / contNbnd.at(m * MRSIZE + n);
      else
        pCorNbnd.at(m * MRSIZE + n) = 0.0;
      COUTPT << pCorNbnd.at(m * MRSIZE + n) << " ";
    }
    COUTPT << endl;
  }
  COUTPT << "\n\noCorr" << endl;
  for(int m = 0; m < MCTOFF; ++m) {
    for(int n = 0; n < MRSIZE; ++n) {
      if(contNbnd.at(m * MRSIZE + n) != 0)
        oCorNbnd.at(m * MRSIZE + n) *= 1.0 / contNbnd.at(m * MRSIZE + n);
      else
        oCorNbnd.at(m * MRSIZE + n) = 0.0;
      COUTPT << oCorNbnd.at(m * MRSIZE + n) << " ";
    }
    COUTPT << endl;
  }
  return 0;
}

void obDelF(world *aSyst, world *cgSys, const int i, const int j,
            const int cgRt, myVect* delF, myCoor& para, myCoor& orth) {
  const partk *atom1, *atom2, *part1, *part2;
  myCoor dist;
  part1 = (*cgSys).obPart(i);
  part2 = (*cgSys).obPart(j);
  myVectZro(TOTM);
  if(!samBnd(part1, part2)) {
    //calculating nonbonded force
    for(int a = 0; a < cgRt; ++a) {
      atom1 = (*aSyst).obPart(i * cgRt + a);
      for(int b = 0; b < cgRt; ++b) {
        atom2 = (*aSyst).obPart(j * cgRt + b);
        (*aSyst).ljonSz(atom1, atom2, dist, UNIM);
        //if(dist < 0.10)
        //  return true;
        (*aSyst).ljonFc(atom1, atom2, 0, dist, UNIM, FORM);
        myVectAdd(TOTM, FORM);
      }
    }
    (*cgSys).ljonSz(part1, part2, dist, UNIM);
    (*cgSys).ljonFc(part1, part2, 2, dist, UNIM, FORM);
    myVectSub(TOTM, FORM);
    myVectCpy(delF, TOTM);
    myVectDot(TOTM, UNIM, &para);
    myVectScl(UNIM, para);
    myVectSub(TOTM, UNIM);
    orth = myVectNrm(TOTM);
  }
  else {
    myVectZro(delF);
    para = 0.0;
    orth = 0.0;
  }
}
