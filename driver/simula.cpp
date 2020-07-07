/****************************************************************************
 * Main driver
 * Default output file is data/trajc.trj
 ***************************************************************************/

#include <iostream>
#include <sstream>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <algorithm>
#include <map>
#include "../libs/const.h"
#include "../libs/ptype.h"
#include "../libs/btype.h"
#include "../libs/atype.h"
#include "../libs/dtype.h"
#include "../libs/rtype.h"
#include "../libs/partk.h"
#include "../libs/bound.h"
#include "../libs/angle.h"
#include "../libs/dihed.h"
#include "../libs/resid.h"
#include "../libs/tlist.h"
#include "../libs/world.h"

using namespace std;
using namespace MD;

tlist types;
world systm;
const vector<ptype> *ptyps;
const vector<btype> *btyps;
const vector<atype> *atyps;
const vector<dtype> *dtyps;
const vector<resid> *rsids;
const vector<partk> *parts;
const vector<bound> *bonds;
const vector<angle> *angls;
const vector<dihed> *dheds;
string space = "   ";
int nptype = 0, nbtype = 0, natype = 0, ndtype = 0, nrtype = 0, 
    npartk = 0, nbound = 0, nangle = 0, ndihed = 0, nresid = 0;

void printVector(const myVect* v) {
  for(int i = 0; i < DIM; ++i)
    cout << myVectGet(v,i) << " ";
}

void printCoor(const vector<partk>& plist) {
  for(int position = 0; position < plist.size(); ++position) {
    cout << "Part number " << position
         << " in the list has position ";
    printVector((plist.at(position)).obPost());
    cout << "and velocity ";
    printVector((plist.at(position)).obVelo());
    cout << endl;
  }
}

void printParameter() {
  ptyps = &(types.obPtypLs());
  btyps = &(types.obBtypLs());
  atyps = &(types.obAtypLs());
  dtyps = &(types.obDtypLs());

  cout << "  ------  " << "In " << paraFl << "  ------  " << endl;
  for(int i = 0; i < (*ptyps).size(); ++i) {
    const ptype *ptyp = &((*ptyps).at(i));
    cout << "  The partk type number " << i << " in the list is "
         << (*ptyp).obName() << " with mass=" << (*ptyp).obMass()
         << ", chrg=" << (*ptyp).obChrg() << ", epsl="
         << (*ptyp).obEpsl() << ", rmin=" << (*ptyp).obRmin()
         << ", sigm=" << (*ptyp).obSigm() << endl;
  }

  for(int i = 0; i < (*btyps).size(); ++i) {
    const btype *btyp = &((*btyps).at(i));
    const ptype *ptyp1 = (*btyp).obLeft();
    const ptype *ptyp2 = (*btyp).obRigh();
    cout << "  The bound type number " << i << " in the list connects "
         << (*ptyp1).obName() << " to " << (*ptyp2).obName()
         << " with relaxing length " << (*btyp).obRelx() << endl;  
  }

  for(int i = 0; i < (*atyps).size(); ++i) {
    const atype *atyp = &((*atyps).at(i));
    const ptype *ptyp1 = (*atyp).obLeft();
    const ptype *ptyp2 = (*atyp).obMidl();
    const ptype *ptyp3 = (*atyp).obRigh();
    cout << "  The angle type number " << i << " in the list connects "
         << (*ptyp2).obName() << " to " << (*ptyp3).obName()
         << " and " << (*ptyp1).obName() << " with relaxing size "
         << ((*atyp).obRelx() / ONEPIE * 180.0) << endl;  
  }

  for(int i = 0; i < (*dtyps).size(); ++i) {
    const dtype *dtyp = &((*dtyps).at(i));
    const ptype *ptyp1 = (*dtyp).obUppr();
    const ptype *ptyp2 = (*dtyp).obLeft();
    const ptype *ptyp3 = (*dtyp).obRigh();
    const ptype *ptyp4 = (*dtyp).obLowr();
    cout << "  The dihed type number " << i << " in the list connects "
         << (*ptyp2).obName() << " and " << (*ptyp3).obName() << " to "
         << (*ptyp1).obName() << " and " << (*ptyp4).obName()
         << " with relaxing size "
         << ((*dtyp).obRelx() / ONEPIE * 180.0) << endl;  
  }
}

//currently bugged
void printStructure() {
  parts = &(systm.obPartLs());
  bonds = &(systm.obBondLs());
  angls = &(systm.obAnglLs());
  dheds = &(systm.obDhedLs());
  rsids = &(systm.obRsidLs());
  int partk_p = 0, bound_p = 0, angle_p = 0, dihed_p = 0;

  cout << "  ------  " << "In " << strcFl << "  ------  " << endl;
  for(int i = 0; i < (*parts).size(); ++i) {
    const partk *atom = &((*parts).at(i));
    const resid *mole = (*atom).obRsid();
    cout << "  The partk number " << i << " in the list has name "
         << (*atom).obName() << " and type "  << (*atom).obMark()
         << " and belongs to resid number " << (*mole).obName()
         << ", with mass "  << (*atom).obMass() << " and charge "
         << (*atom).obChrg() << endl;  
  }

  for(int i = 0; i < (*bonds).size(); ++i) {
    const bound *bond = &((*bonds).at(i));
    const partk *atom1 = (*bond).obLeft();
    const partk *atom2 = (*bond).obRigh();
    const resid *mole = (*bond).obRsid();
    cout << "  The bound number " << i << " in the list connects "
         << (*atom1).obName() << " to " << (*atom2).obName()
         << " - which belong to resid number " << (*mole).obName()
         << " - with relaxing length " << (*bond).obRelx() << endl;  
  }

  for(int i = 0; i < (*angls).size(); ++i) {
    const angle *angl = &((*angls).at(i));
    const partk *atom1 = (*angl).obLeft();
    const partk *atom2 = (*angl).obMidl();
    const partk *atom3 = (*angl).obRigh();
    const resid *mole = (*angl).obRsid();
    cout << "  The angle number " << i+1 << " in the list connects "
         << (*atom2).obName() << " to " << (*atom3).obName()
         << " and " << (*atom1).obName()
         << " - which belong to resid number " << (*mole).obName()
         << " - with relaxing size " << (*angl).obRelx() << endl;  
  }

  for(int i = 0; i < (*dheds).size(); ++i) {
    const dihed *dhed = &((*dheds).at(i));
    const partk *atom1 = (*dhed).obUppr();
    const partk *atom2 = (*dhed).obLeft();
    const partk *atom3 = (*dhed).obRigh();
    const partk *atom4 = (*dhed).obLowr();
    const resid *mole = (*dhed).obRsid();
    cout << endl;
  }

  for(int i = 0; i < (*rsids).size(); ++i) {
    const resid *mole = &((*rsids).at(i));
    cout << "  The resid number " << i << " in the list has type "
         << (*((*mole).obRtyp())).obName() << ", containing:"
         << endl << space << (*mole).obNpartk() << " parts:"
         << endl;  
    for(int j = 0; j < (*mole).obNpartk(); ++j) {
      const partk *atom = (*mole).obPart(j+1);
      cout << space << "  The partk number " << j << " has type "
           << (*((*mole).obPart(j+1))).obMark() << " and name "
           << ((*parts).at(j+partk_p)).obName() << endl;  
    }
    cout << endl << space << (*mole).obNbound() << " bonds:" << endl;
    for(int j = 0; j < (*mole).obNbound(); ++j) {
      const bound *bond = (*mole).obBond(j+1);
      cout << space << "  The bound number " << j << " connects "
           << (*((*bond).obLeft())).obName() << " to "
           << (*((*bond).obRigh())).obName()
           << ", with relaxing length "
           << ((*bonds).at(j+bound_p)).obRelx() << endl;  
    }
    cout << endl << space << (*mole).obNangle() << " angls:"
         << endl; 
    cout << endl << space << (*mole).obNdihed() << " dheds:"
         << endl; 
    partk_p += (*mole).obNpartk();
    bound_p += (*mole).obNbound();
  }    
} 

void printCoordinate() {
  parts = &(systm.obPartLs());

  cout << "  ------  " << "In " << coorFl << "  ------  " << endl;

  for(int i = 0; i <(*parts).size(); ++i) {
    cout << "  The partk number " << i << " in the list has coor ";
    printVector(((*parts).at(i)).obPost());
    cout << " and velocity ";
    printVector(((*parts).at(i)).obVelo());
    cout << endl;
  }
}

int main() {
  int duration, relxtion, interval;
  char optTrj[256];
  myCoor sTemp, sDcon, sIcon, ener;
  myVect  *bVec, *v1, *v2;

  bVec = myVectMalloc(DIM);
  v1 = myVectMalloc(DIM);
  v2 = myVectMalloc(DIM);
  parts = &(systm.obPartLs());

  //Initialize simulation conditions
  systm.mkCond(COUTPT, CINPUT);

  COUTPT << "Choose the file path (leaving blank gives default path)" << endl;
  CINPUT.ignore();
  CINPUT.getline(optTrj, 256);
  if(strlen(optTrj) == 0)
    strcpy(optTrj, trajFl);
  COUTPT << optTrj << endl;

  COUTPT << "Choose the number of steps for the relaxation" << endl;
  CINPUT >> relxtion;
  COUTPT << relxtion << endl;

  COUTPT << "Choose the number of steps for the simulation" << endl
       << "One step is " << STPSZE << " picoseconds" << endl;
  CINPUT >> duration;
  COUTPT << duration << endl;

  COUTPT << "Choose the interval at which the trajectory is recorded" << endl;
  CINPUT >> interval;
  COUTPT << interval << endl;

  //Read data
  COUTPT << "Reading parameter file..." << endl;
  ifstream para(paraFl, ios_base::in);
  types.adType(para);
  para.close();
  //printParameter();

  COUTPT << "Reading topology file..." << endl;
  ifstream topo(topoFl, ios_base::in);
  types.adType(topo);
  topo.close();

  COUTPT << "Reading structure file..." << endl;
  ifstream strc(strcFl, ios_base::in);
  systm.iniSys(strc, &types);
  strc.close();
  //printStructure();

  ifstream coor(coorFl, ios_base::in);
  systm.iniTrj(coor);
  coor.close();
  //printCoordinate();

  systm.adMemo();
  ifstream memory(memoFl, ios_base::in);
  if(systm.obDcon() >= 5 && systm.obDcon() < 11)
    systm.iniMem(memory);

  //Start simulation
  ofstream trajectory(optTrj, ios_base::out);
  systm.wrCond(trajectory);
  trajectory << duration << " " << interval << "\n";
  srand(time(NULL));
  cout << "Seed: " << ENGINE << endl;
  if(systm.obIcon() == 32) {
    cout << "Thermalizing the system..." << endl;
    systm.mkIcon(31);
    for(int i = 0; i < relxtion; ++i) {
      systm.nwTraj();
      if(i % interval == 0) { 
        cout << optTrj << " -> " << i * STPSZE << " ";
        ener = systm.obKine();
        cout << "KE=" << ener << " "
             << "T=" << ener / (systm.obNpartk() * DIM / 2) << endl; 
      }
    }
    cout << "Swtiching back to self-consistent DPD-VV." << endl
         << "Temperature is " << systm.obSysT() << endl;
    systm.mkIcon(32);
  }
  /*
  if(systm.obDcon() >= 8 && systm.obDcon() < 11) {
    cout << "Thermalizing the system..." << endl;
    sDcon = systm.obDcon();
    systm.mkDcon(2);
    for(int i = 0; i < relxtion; ++i) {
      systm.nwTraj();
      if(i % interval == 0) { 
        cout << optTrj << " -> " << i * STPSZE << " ";
        ener = systm.obKine();
        cout << "KE=" << ener << " "
             << "T=" << ener / (systm.obNpartk() * DIM / 2) << endl; 
      }
    }
    cout << "Swtiching back to auxiliary MZ." << endl
         << "Temperature is " << systm.obSysT() << endl;
    systm.mkDcon(sDcon);
  }
  */
  for(int i = 0; i < (relxtion + duration); ++i) {
    if(i % interval == 0) {
      if(i >= relxtion)
        systm.wrTraj(trajectory);
      cout << optTrj << " -> t=" << systm.obTime() << " ";
      ener = systm.obKine();
      cout << "KE=" << ener << " "
          << "T/<T>=" << ener / (systm.obNpartk() * systm.obTemp() * BOLTZM
                              * DIM / 2) << " "
           << "p=" << myVectNrm(systm.obMomt()) << " ";
      cout << endl;
    }
    systm.nwTraj();
  }
  trajectory.close();
  return 0;
} 
