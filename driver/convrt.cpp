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
#include "../libs/tlist.h"
#include "../libs/world.h"

using namespace MD;
using namespace std;

int obPtypNo(const ptype* ptyp, const tlist* typs);
int obBtypNo(const btype* btyp, const tlist* typs);
int obAtypNo(const atype* atyp, const tlist* typs);
int obDtypNo(const dtype* dtyp, const tlist* typs);

int main(int argc, char *argv[]) {
  world systm;
  tlist types;
  const vector<partk> *parts;
  char optFle[256], outFle[256], inpFle[256], format[256];
  int ctrl;
  FILE *pf;

  cout << "This program converts the output from simula.exe to standard\n"
       << "output. Currently support LAMMPS dump format, GROMACS top and gro\n"
       << "format\n"
       << endl;

  cout << "Choose the format for conversion:\n"
       << "1 for conversion from .gro, 2 for conversion to .gro,\n"
       << "3 for conversion from .top, 4 for conversion to .top,\n"
       << "5 for conversion from .dump, 6 for conversion to .dump\n"
       << endl;
  cin >> ctrl;

  if(argc != 3) {
    cerr << "Usage:" << argv[0] << " <inpt> <outp>" << endl;
    cerr << "Using default path" << endl;
    if(ctrl == 3) {
      strcpy(inpFle, topoFl);
      strcat(inpFle, ".top");
      strcpy(outFle, topoFl);
    }
    else if(ctrl == 6) {
      strcpy(inpFle, trajFl);
      strcpy(outFle, trajFl);
      strcat(outFle, ".dump");
    }
  }
  else {
    strcpy(inpFle, argv[1]);
    strcpy(outFle, argv[2]);
  }
  cout << argv[0] << " < " << inpFle << " > " << outFle << endl;

  ifstream para(paraFl, ios_base::in);
  types.adType(para);
  para.close();
  ifstream topo(topoFl, ios_base::in);
  types.adType(topo);
  topo.close();
  ifstream strc(strcFl, ios_base::in);
  systm.iniSys(strc, &types);
  strc.close();   
  ifstream coor(coorFl, ios_base::in);
  systm.iniTrj(coor);
  coor.close();

  parts = &(systm.obPartLs());

  if(ctrl == 1) {
  }
  else if(ctrl == 3) {
    string dummy;
    vector <string> rName, pName;
    vector <int> pNums, bNums, aNums, dNums;
    ifstream inpt(inpFle, ios_base::in);
    ifstream topo(topoFl, ios_base::in);
    ofstream dump("scratch.txt", ios_base::out);

    while(inpt >> dummy) {
      if(dummy == ";")
        inpt.ignore(256, '\n');
      else if(dummy == "moleculetype") {
        inpt.ignore(256, '\n');
        inpt.ignore(256, '\n');
        inpt >> dummy;
        rName.push_back(dummy);
        
        pNums.push_back(0);
        bNums.push_back(0);
        aNums.push_back(0);
        dNums.push_back(0);
        inpt.ignore(256, '\n');
      }
    }
    while(!topo.eof()) {
      topo.getline(format, 256);
      if(strcmp(format, "END"))
        dump << format << endl;
      else
        break;
    }

  }
  else if(ctrl == 6) {
    int duration, interval;
    //Write structure/topology
    ofstream file;
    file.open(string(defNme) + "struc.lj");
    file << "LAMMPS input file" << endl;
    file << endl;
    file << systm.obNpartk() << " atoms" << endl;
    file << systm.obNbound() << " bonds" << endl;
    file << systm.obNangle() << " angles" << endl;
    file << systm.obNdihed() << " dihedrals" << endl;
    file << "0 impropers" << endl;
    file << endl;
    file << types.obNptype() << " atom types" << endl;
    if(types.obNbtype() != 0)
      file << types.obNbtype() << " bond types" << endl;
    if(types.obNatype() != 0)
      file << types.obNatype() << " angle types" << endl;
    if(types.obNdtype() != 0)
      file << types.obNdtype() << " dihedral types" << endl;
    file << endl;
    file << "0.000 " << myVectGet(systm.obEdge(), 0)
         << " xlo xhi" << endl;
    file << "0.000 " << myVectGet(systm.obEdge(), 1)
         << " ylo yhi" << endl;
    file << "0.000 " << myVectGet(systm.obEdge(), 2)
         << " zlo zhi" << endl;
    file << endl;
    file << "Masses" << endl;
    file << endl;
    for(int i = 0; i < types.obNptype(); ++i)
      file << (i + 1) << " " << (*(types.obPtyp(i))).obMass() << endl;
    file << endl;
    file << "Nonbond Coeffs" << endl;
    file << endl;
    for(int i = 0; i < types.obNptype(); ++i)
      file << (i + 1) << " " << (*(types.obPtyp(i))).obEpsl() << " "
           << (*(types.obPtyp(i))).obRmin() << endl;
    file << endl;
    file << "Bond Coeffs" << endl;
    file << endl;
    for(int i = 0; i < types.obNbtype(); ++i)
      file << (i + 1) << " " << (*(types.obBtyp(i))).obRelx() << " "
           << (*(types.obBtyp(i))).obRigd() << endl;
    file << endl;
    file << "Angle Coeffs" << endl;
    file << endl;
    for(int i = 0; i < types.obNatype(); ++i)
      file << (i + 1) << " " << (*(types.obAtyp(i))).obRelx() << " "
           << (*(types.obAtyp(i))).obRigd() << endl;
    file << endl;
    file << "Dihedral Coeffs" << endl;
    file << endl;
    for(int i = 0; i < types.obNdtype(); ++i)
      file << (i + 1) << " " << (*(types.obDtyp(i))).obRelx() << " "
           << (*(types.obDtyp(i))).obRig1() << " "
           << (*(types.obDtyp(i))).obRig2() << " "
           << (*(types.obDtyp(i))).obRig3() << endl;
    file << endl;
    file << "Atoms" << endl;
    file << endl;
    for(int i = 0; i < systm.obNpartk(); ++i) {
      file << (i + 1) << " " 
           << (*((*(systm.obPart(i))).obRsid())).obName() << " "
           << obPtypNo((*(systm.obPart(i))).obType(), &types) << " "
           //<< (*(systm.obPart(i))).obChrg()
           << " ";
      for(int k = 0; k < DIM; ++k)
        file << myVectGet((*(systm.obPart(i))).obPost(), k) << " ";
      file << endl;
    }
    file << endl;
    file << "Velocities" << endl;
    file << endl;
    for(int i = 0; i < systm.obNpartk(); ++i) {
      file << (i + 1) << " "; 
      for(int k = 0; k < DIM; ++k)
        file << myVectGet((*(systm.obPart(i))).obVelo(), k) << " ";
      file << endl;
    }
    file << endl;
    file << "Bonds" << endl;
    for(int i = 0; i < systm.obNbound(); ++i) {
      file << (i + 1) << " " 
           << (obBtypNo((*(systm.obBond(i))).obType(), &types) + 1) << " "
           << (*((*(systm.obBond(i))).obLeft())).obName() << " "
           << (*((*(systm.obBond(i))).obRigh())).obName() << endl;
    }
    file << endl;
    file << "Angle" << endl;
    file << endl;
    for(int i = 0; i < systm.obNangle(); ++i) {
      file << (i + 1) << " " 
           << (obAtypNo((*(systm.obAngl(i))).obType(), &types) + 1) << " "
           << (*((*(systm.obAngl(i))).obLeft())).obName() << " "
           << (*((*(systm.obAngl(i))).obMidl())).obName() << " "
           << (*((*(systm.obAngl(i))).obRigh())).obName() << endl;
    }
    file << endl;
    file << "Dihedrals" << endl;
    file << endl;
    for(int i = 0; i < systm.obNdihed(); ++i) {
      file << (i + 1) << " " 
           << (obDtypNo((*(systm.obDhed(i))).obType(), &types)  + 1) << " "
           << (*((*(systm.obDhed(i))).obUppr())).obName() << " "
           << (*((*(systm.obDhed(i))).obLeft())).obName() << " "
           << (*((*(systm.obDhed(i))).obRigh())).obName() << " "
           << (*((*(systm.obDhed(i))).obLowr())).obName() << endl;
    }
    file << endl;
    file << "Impropers" << endl;
    file << endl;
    file.close();
    //Write trajectory
    ifstream traj(inpFle, ios_base::in);
    pf = fopen(outFle, "w");
    systm.iniCon(traj);
    systm.mkDcon(0);
    traj >> duration >> interval;
    for(int i = 0; i < (duration / interval); ++i) {
      systm.iniTrj(traj);
      fprintf(pf, "ITEM: TIMESTEP\n");
      fprintf(pf, crdFmt, systm.obTime());
      fprintf(pf, "\n");
      fprintf(pf, "ITEM: NUMBER OF ATOMS\n");
      fprintf(pf, "%-80d\n", systm.obNpartk());
      fprintf(pf, "ITEM: BOX BOUNDS pp pp pp\n");
      switch(systm.obBcon()) {
        case 1:
          strcpy(format, crdFmt);
          strcat(format, crdFmt);
          strcat(format, "\n");
          for(int k = 0; k < DIM; ++k)
            fprintf(pf, format,
                        0.0, myVectGet(systm.obEdge(), k));
          break;
      }
      fprintf(pf, "ITEM: ATOMS id mol type x y z vx vy vz\n");
      for(int j = 0; j < systm.obNpartk(); ++j) {
        strcpy(format, sNmFmt);
        strcat(format, sNmFmt);
        strcat(format, sNmFmt);
        fprintf(pf, format,
                    (((*parts).at(j)).obName()).c_str(),
                    ((*(((*parts).at(j)).obRsid())).obName()).c_str(),
                    (((*parts).at(j)).obMark()).c_str());
        for(int k = 0; k < DIM; ++k)
          fprintf(pf, crdFmt, myVectGet(((*parts).at(j)).obPost(), k));
        for(int k = 0; k < DIM; ++k)
          fprintf(pf, crdFmt, myVectGet(((*parts).at(j)).obVelo(), k));
        fprintf(pf, "\n");
      }
    }
  }
  return 0;
}

int obPtypNo(const ptype* ptyp, const tlist* typs) {
  const vector<ptype>* ptyps = &((*typs).obPtypLs());
  for(int i = 0; i < (*typs).obNptype(); i++)
    if((*ptyps).at(i) == (*ptyp))
      return i;
  cerr << "Cannot find ptype" << endl;
}

int obBtypNo(const btype* btyp, const tlist* typs) {
  const vector<btype>* btyps = &((*typs).obBtypLs());
  for(int i = 0; i < (*typs).obNbtype(); i++)
    if((*btyps).at(i) == (*btyp))
      return i;
  cerr << "Cannot find btype" << endl;
}

int obAtypNo(const atype* atyp, const tlist* typs) {
  const vector<atype>* atyps = &((*typs).obAtypLs());
  for(int i = 0; i < (*typs).obNatype(); i++)
    if((*atyps).at(i) == (*atyp))
      return i;
  cerr << "Cannot find atype" << endl;
}

int obDtypNo(const dtype* dtyp, const tlist* typs) {
  const vector<dtype>* dtyps = &((*typs).obDtypLs());
  for(int i = 0; i < (*typs).obNdtype(); i++)
    if((*dtyps).at(i) == (*dtyp))
      return i;
  cerr << "Cannot find dtype" << endl;
}
