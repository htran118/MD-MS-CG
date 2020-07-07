/****************************************************************************
 * Coarse-graining a system
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

world systm;
tlist types;
vector<resid> *rsids;
vector<partk> *parts;
vector<bound> *bonds;
vector<angle> *angls;
vector<dihed> *dheds;

bool isSequ(vector <int> *intLs) {
  if((*intLs).size > 0) {
    sort((*intLs).begin(), (*intLs).end());
    for(int i = 0; i < (*intLs).size() - 1; ++i)
      if((*intLs).at(i + 1) - (*intLs).at(i) != 1) {
        cerr << "The input is not sequential" << endl;
        return false;
      }
    return true;
  }
  else
    return false;
}

void makeCG(const string rName) {
  vector <resid*> oRsid;
  vector <resid*> nRsid;
  vector <int> CGlst;
  int ctrl, dumy;
  string inps, name;
  bool same;

  for(int i = 0; i < (*rsids).size(); ++i)
    if(((*rsids).at(i)).obRtyp() == types.obRtyp(rName))
      oRsid.push_back(const_cast <resid*> (&((*rsids).at(i))));
  mole = oRsid.at(0);
  cout << "Residue type " << rName << " contains "
       << (*mole).obNpartk() << " atoms:\n";
  for(int i = 0; i < (*mole).obNpartk(); ++i)
    cout << "Partk " << i << " has mark "
         << (*((*mole).obPart(i))).obMark() << endl;
  cout << "Choose the number of atoms to coarse-grain\n";
  cin >> ctrl;
  cout << "Choose the atoms to coarse-grain\n";
  if(isSequ(&CGlst) {
    CGlst.clear();
    for(int i = 0; i < ctrl; ++i)
      cin >> dumy;
      CGlst.push_back(dumy);
    }
  }
  strcpy(
    
  






}

int main(int argc, char *argv[]) {
  vector <string> rName;
  vector <int> rNums;
  char optFle[256], outFle[256], inpFle[256], format[256];
  int ctrl;
  char inpt;
  FILE *pf;
  string inps;
  bool same, again;

  cout << "This program create a CG model from a given system\n";
  if(argc != 3) {
    cerr << "Usage:" << argv[0] << " <inpt> <outp>" << endl;
    cerr << "Using default path" << endl;
    strcpy(inpFle, coorFl);
    strcpy(outFle, coorFl);
    strcat(outFle, "CG");
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
  rsids = const_cast <vector <resid>> (&(systm.obRsidLs()));
  parts = const_cast <vector <partk>> (&(systm.obPartLs()));
  bonds = const_cast <vector <bound>> (&(systm.obBondLs()));
  angls = const_cast <vector <angle>> (&(systm.obAnglLs()));
  dheds = const_cast <vector <dihed>> (&(systm.obDhedLs()));

  for(int i = 0; i < (*rsids).size(); ++i) {
    inps = (*(((*rsids).at(i)).obRtyp())).obName();
    same = false;
    for(int j = 0; j < rName.size() && !same; ++j)
      if(same = (inps == rName.at(j)))
        rNums.at(j) += 1;
    if(!same) {
      rName.push_back(inps);
      rNums.push_back(1);
    }
  }
  cout << "There are " << (*rsids).size() << " residues, of which:\n";
  for(int i = 0; i < rName.size(); ++i)
    cout << rNums.at(i) << " are of rtype " << rName.at(i) << endl;
   
  again = false;
  while(!again) {
    cout << "Choose the rtype to coarse-grain\n";
    cin >> inps;
    same = false;
    while(!same) {
      for(int i = 0; i < rName.size() && !same; ++i)
        if(same = (inps == rName.at(i)))
          makeCG(inps);
      if(!same) {
        cout << "The given rtype " << inps << " is not found.\n";
        cout << "Choose the rtype to coarse-grain\n";
        cin >> inps;
      }
    }
    cout << "Coarse-grain another type?\n";
    cin >> inpt;
    again = (inpt == 'y');
  }
  return 0;
}
