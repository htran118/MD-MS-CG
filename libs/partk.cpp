#include <cassert>
#include <algorithm>
#include "partk.h"

using namespace std;

namespace MD {
  //CONSTRUCTORS
  partk::partk() {
    post = myVectMalloc(DIM);
    velo = myVectMalloc(DIM);
    accl = myVectMalloc(DIM);
    name = mark = emptyNme;
  }

  partk::partk(const string sMark, const ptype* sType) {
    post = myVectMalloc(DIM);
    velo = myVectMalloc(DIM);
    accl = myVectMalloc(DIM);
    mark = sMark;
    type = sType;
  }

  partk::partk(const partk& sAtom) {
    post = myVectMalloc(DIM);
    velo = myVectMalloc(DIM);
    accl = myVectMalloc(DIM);
    myVectCpy(post, sAtom.obPost());
    myVectCpy(velo, sAtom.obVelo());
    myVectCpy(accl, sAtom.obAccl());
    mole = sAtom.obRsid();
    name = sAtom.obName();
    mark = sAtom.obMark();
    type = sAtom.obType();
  }

  //BUG: input reading does not support changing partk type
  partk::partk(istream& iput, const world* syst) {
    string rNme;
    iput >> name >> rNme >> mark;
    post = myVectMalloc(DIM);
    velo = myVectMalloc(DIM);
    accl = myVectMalloc(DIM);
    mole = (*syst).obRsid(rNme);
    if(mole != NULL) 
      type = (*((*mole).obRtyp())).obPtyp(mark);
    else
      cerr << "Partk " << mark << " does not belong to any resid."
           << endl;
  }

  partk::~partk() {
    myVectMemdel(post);
    myVectMemdel(velo);
    myVectMemdel(accl);
  }

  //MODIFIERS
  const partk& partk::operator =(const partk& sAtom) {
    if(this != &sAtom) {
      myVectCpy(post, sAtom.obPost());
      myVectCpy(velo, sAtom.obVelo());
      myVectCpy(accl, sAtom.obAccl());
      mark = sAtom.obMark();
      type = sAtom.obType();
    }
    return *this;
  }

  void partk::mkForc(const myVect* sForc) {
    myVectCpy(accl, sForc);
    myVectInv(accl, (*type).obMass());
  }

  //BUG: doesn't seem to have a good way to implement
  void partk::adForc(const myVect* fDif) {
    cerr << "Not supporting directing adding force to partk" << endl;
  }

  //add bound
  void partk::adBond(bound* bond) {
    //NOTE: physically should be vale, but some parts contain more bonds
    //than possible to accomodate advanced techniques
    assert(bonds.size() < MAXNBOND);
    bonds.push_back(bond);
  }

  //remove bound
  void partk::rmBond(bound* bond) {
    bonds.erase(remove(bonds.begin(), bonds.end(), bond), bonds.end());
    //NOTE: maybe I can do this?
    //vector<bound*>::iterator pos = 
    //find(bonds.begin(), bonds.end(), bond);
    //if(pos != bonds.end())
    //  bonds.erase(pos);
  }

  void partk::rmBond(const int indx) {
    assert(indx >= 0 && indx < bonds.size());
    bonds.erase(bonds.begin() + indx);
  }

  //add angle
  void partk::adAngl(angle* angl) {
    assert(angls.size() < MAXNANGL);
    angls.push_back(angl);
  }

  //remove angle
  void partk::rmAngl(angle* angl) {
    angls.erase(remove(angls.begin(), angls.end(), angl), angls.end());
    //NOTE: maybe I can do this?
    //vector<angle*>::iterator pos = 
    //find(angls.begin(), angls.end(), angl);
    //if(pos != angls.end())
    //  angls.erase(pos);
  }

  void partk::rmAngl(const int indx) {
    assert(indx >= 0 && indx < angls.size());
    angls.erase(angls.begin() + indx);
  }

  //add dihed
  void partk::adDhed(dihed* dhed) {
    assert(dheds.size() < MAXNANGL);
    dheds.push_back(dhed);
  }

  //remove dihed
  void partk::rmDhed(dihed* dhed) {
    dheds.erase(remove(dheds.begin(), dheds.end(), dhed), dheds.end());
    //NOTE: maybe I can do this?
    //vector<dihed*>::iterator pos = 
    //find(dheds.begin(), dheds.end(), dhed);
    //if(pos != dheds.end())
    //  dheds.erase(pos);
  }

  void partk::rmDhed(const int indx) {
    assert(indx >= 0 && indx < dheds.size());
    dheds.erase(dheds.begin() + indx);
  }
 
  //initialize trajectory
  void partk::iniTrj(istream& iput) {
    myCoor n;
    for(int i = 0; i < DIM; ++i) {
      iput >> n;
      myVectSet(post, i, n);
    }
    for(int i = 0; i < DIM; ++i) {
      iput >> n;
      myVectSet(velo, i, n);
    }
  }

  //update using Verlet integrator
  void partk::verlet() {
    //original Verlet
    //velo is x(t-dt)
    //x(t+dt) = 2x(t) - x(t-dt) + a(t)dt2
    myVectScl(accl, STEPSQ);
    myVectSub(accl, velo);
    myVectCpy(velo, post);
    myVectScl(post, 2);
    myVectAdd(post, accl);
  }

  //update using leap frog integrator
  void partk::lpFrgA() {
    //leapfrog
    //velo is v(t-0.5dt)
    //v(t+0.5dt) = v(t-0.5dt) + a(t)dt
    myVectScl(accl, STPSZE);
    myVectAdd(velo, accl);
    //xprime = x(t) + 0.5*v(t+0.5dt)dt
    myVectCpy(accl, velo);
    myVectScl(accl, STPHLF);
    myVectAdd(post, accl);
  }

  void partk::lpFrgB() {
    //leapfrog for DPD
    //x(t+dt) = xprime + 0.5*vprime*dt
    myVectCpy(accl, velo);
    myVectScl(accl, STPHLF);
    myVectAdd(post, accl);
  }

  //update using velocity Verlet integrator
  void partk::vVerlA() {
    //velocity Verlet 1
    //v(t+0.5dt) = v(t) + 0.5a(t)dt
    myVectScl(accl, STPHLF);
    myVectAdd(velo, accl);
    //x(t+dt) = x(t) + v(t+0.5dt)dt
    myVectCpy(accl, velo);
    myVectScl(accl, STPSZE);
    myVectAdd(post, accl);
  }

  void partk::vVerlB() {
    //velocity Verlet 2
    //v(t) = v(t+0.5dt) + 0.5a(t+dt)dt
    myVectScl(accl, STPHLF);
    myVectAdd(velo, accl);
    myVectScl(accl, STPIHF);
  }

  //reflect velocity (for hard-sphere collision)
  void partk::rflect(const int indx) {
    myVectSet(velo, indx, -myVectGet(velo, indx));
  }

  //OBSERVERS
  //write trajectory
  void partk::wrTraj(ostream& oput) const {
    //write out update trajectories
    for(int i = 0; i < DIM; ++i)
      oput << myVectGet(post, i) << '\t';
    for(int i = 0; i < DIM; ++i)
      oput << myVectGet(velo, i) << '\t';
    oput << endl;
  }

  //write position
  void partk::wrCoor(ostream& oput) const {
    //write out update coordinates
    for(int i = 0; i < DIM; ++i)
      oput << myVectGet(post, i) << '\t';
    oput << endl;
  }

  //obtain position
  myCoor partk::obCoor(const int indx) const {
    assert(indx < DIM);
    return myVectGet(post, indx);
  }

  //obtain bond
  const bound* partk::obBond(const int indx) const {
    assert(indx >= 0 && indx < bonds.size());
    return bonds.at(indx);
  }

  //obtain number of bond
  int partk::obBondNo(const bound* bond) const {
    int pos = find(bonds.begin(), bonds.end(), bond) - bonds.begin();
    if(pos < bonds.size())
      return (pos + 1);
    else
      return MAXNBOND;
  }

  //obtain angle
  const angle* partk::obAngl(const int indx) const {
    assert(indx >= 0 && indx < angls.size());
    return angls.at(indx);
  }

  //obtain number of angle
  int partk::obAnglNo(const angle* angl) const {
    int pos = find(angls.begin(), angls.end(), angl) - angls.begin();
    if(pos < angls.size())
      return (pos + 1);
    else
      return MAXNANGL;
  }

  //obtain dihed
  const dihed* partk::obDhed(const int indx) const {
    assert(indx >= 0 && indx < dheds.size());
    return dheds.at(indx);
  }

  //obtain number of dihed
  int partk::obDhedNo(const dihed* dhed) const {
    int pos = find(dheds.begin(), dheds.end(), dhed) - dheds.begin();
    if(pos < dheds.size())
      return (pos + 1);
    else
      return MAXNANGL;
  }

  //NON-MEMBER FUNCTIONS
  //obtain position difference
  //the result vector points from atom1 to atom2
  void obXdif(const partk* atom1, const partk* atom2,
              myVect* const xDif)  {
    myVectCpy(xDif, (*atom2).obPost());
    myVectSub(xDif, (*atom1).obPost());
  }

  //obtain velocity difference
  void obVdif(const partk* atom1, const partk* atom2,
              myVect* const vDif)  {
    myVectCpy(vDif, (*atom2).obVelo());
    myVectSub(vDif, (*atom1).obVelo());
  }

  const bool samRsd(const partk* atom1, const partk* atom2) {
    return ((*atom1).obRsid() == (*atom2).obRsid());
  }

  //check if 2 parts are connected directly through any topology structure
  //such as bond, angle, dihedral angle
  //used for calculating non-bonded interaction
  const bool samTop(const partk* atom1, const partk* atom2) {
    bool isBond = false;
    if(BONDFORC > 2 && !isBond)
      isBond |= samDih(atom1, atom2);
    if(BONDFORC > 1 && !isBond)
      isBond |= samAng(atom1, atom2);
    if(BONDFORC > 0 && !isBond)
      isBond |= samBnd(atom1, atom2);
    return isBond;
  }

  const bool samDih(const partk* atom1, const partk* atom2) {
    bool isBond = false;
    const dihed* dhed;
    for(int i = 0; i < (*atom1).obNdihed(); ++i) {
      dhed = (*atom1).obDhed(i);
      isBond |= ((*dhed).obUppr() == atom2 || (*dhed).obLeft() == atom2 ||
                 (*dhed).obRigh() == atom2 || (*dhed).obLowr() == atom2);
    }
    return isBond;
  }

  const bool samAng(const partk* atom1, const partk* atom2) {
    bool isBond = false;
    const angle* angl;
    for(int i = 0; i < (*atom1).obNangle(); ++i) {
      angl = (*atom1).obAngl(i);
      isBond |= ((*angl).obLeft() == atom2 || (*angl).obMidl() == atom2 ||
                 (*angl).obRigh() == atom2);
    }
    return isBond;
  }

  const bool samBnd(const partk* atom1, const partk* atom2) {
    bool isBond = false;
    for(int i = 0; i < (*atom1).obNbound(); ++i)
      isBond |= ((*((*atom1).obBond(i))).obLeft() == atom2 ||
                 (*((*atom1).obBond(i))).obRigh() == atom2);
    return isBond;
  }

  //obtain bound
  const bound* obBond(const partk* atom1, const partk* atom2) {
    const bound* bond;
    for(int i = 0; i < (*atom1).obNbound(); ++i) {
      bond = (*atom1).obBond(i);
      if((*bond).obLeft() == atom2 || (*bond).obRigh() == atom2)
        return bond;
    }
    cerr << "Cannot find bond connecting " << (*atom1).obName()
         << " to " << (*atom2).obName() << endl;
    return NULL;
  }

  //obtain angle
  //BUG: currently only return the correct angle if the parts are given
  //in correct order
  const angle* obAngl(const partk* atom1, const partk* atom2,
                      const partk* atom3) {
    const angle* angl;
    for(int i = 0; i < (*atom2).obNangle(); ++i) {
      cout << (*atom2).obNangle() << endl;
      angl = (*atom2).obAngl(i);
      if(((*angl).obLeft() == atom1 && (*angl).obMidl() == atom2 &&
          (*angl).obRigh() == atom3) ||
         ((*angl).obLeft() == atom3 && (*angl).obMidl() == atom2 &&
          (*angl).obRigh() == atom1))
        return angl;
    }
    cerr << "Cannot find angle connecting " << (*atom2).obName()
         << " to " << (*atom1).obName()
         << " and " << (*atom3).obName() << endl;
    return NULL;
  }

  //obtain dihed
  //BUG: currently only return the correct angle if the parts are given
  //in correct order
  const dihed* obDhed(const partk* atom1, const partk* atom2,
                      const partk* atom3, const partk* atom4) {
    const dihed* dhed;
    for(int i = 0; i < (*atom1).obNdihed(); ++i) {
      dhed = (*atom2).obDhed(i);
      if(((*dhed).obUppr() == atom1 && (*dhed).obLeft() == atom2 &&
          (*dhed).obRigh() == atom3 && (*dhed).obLowr() == atom4) ||
         ((*dhed).obUppr() == atom4 && (*dhed).obLeft() == atom3 &&
          (*dhed).obRigh() == atom2 && (*dhed).obLowr() == atom1))
        return dhed;
    }
    cerr << "Cannot find dihedral angle connecting " << (*atom2).obName()
         << " and " << (*atom3).obName()
         << " to " << (*atom1).obName()
         << " and " << (*atom4).obName() << endl;
    return NULL;
  }

}
