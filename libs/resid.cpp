#include <cassert>
#include <algorithm>
#include "resid.h"
using namespace std;

namespace MD {
  //CONSTRUCTORS
  resid::resid(istream& iput, const tlist* types) {
    string rName;
    iput >> name >> rName;
    type = (*types).obRtyp(rName);
  }

  //MODIFIERS
  void resid::adPart(partk* atom) {
    assert(parts.size() < MAXNPART);
    parts.push_back(atom);
  }

  void resid::adBond(bound* bond) {
    assert(bonds.size() < MAXNBOND);
    bonds.push_back(bond);
  }

  void resid::adAngl(angle* angl) {
    assert(angls.size() < MAXNANGL);
    angls.push_back(angl);
  }

  void resid::adDhed(dihed* dhed) {
    assert(dheds.size() < MAXNDHED);
    dheds.push_back(dhed);
  }

  void resid::rmPart(partk* atom) {
    parts.erase(remove(parts.begin(), parts.end(), atom), parts.end());
  }

  void resid::rmBond(bound* bond) {
    bonds.erase(remove(bonds.begin(), bonds.end(), bond), bonds.end());
  }

  void resid::rmAngl(angle* angl) {
    angls.erase(remove(angls.begin(), angls.end(), angl), angls.end());
  }

  void resid::rmDhed(dihed* dhed) {
    dheds.erase(remove(dheds.begin(), dheds.end(), dhed), dheds.end());
  }

  void resid::rmBond(const int indx) {
    assert(indx >= 0 && indx < parts.size());
    parts.erase(parts.begin() + indx);
  }

  void resid::rmPart(const int indx) {
    assert(indx >= 0 && indx < bonds.size());
    bonds.erase(bonds.begin() + indx);
  }

  void resid::rmAngl(const int indx) {
    assert(indx >= 0 && indx < angls.size());
    angls.erase(angls.begin() + indx);
  }

  void resid::rmDhed(const int indx) {
    assert(indx >= 0 && indx < dheds.size());
    dheds.erase(dheds.begin() + indx);
  }

  const resid& resid::operator =(const resid& mole) {
    if(this != &mole) {
      name = mole.obName();
      //assign the components of the residue
    }
    return *this;
  }

  void resid::rflect(const int indx) {
    for(int i = 0; i < parts.size(); ++i)
      (*(parts.at(i))).rflect(indx);
  }

  //OBSERVERS
  const partk* resid::obPart(const int indx) const {
    assert(indx >= 0 && indx < parts.size());
    return parts.at(indx);
  }

  int resid::obPartNo(partk* atom) const {
    int post = find(parts.begin(), parts.end(), atom) - parts.begin();
    if(post < parts.size())
      return (post + 1);
    else
      return MAXNPART;
  }

  const bound* resid::obBond(const int indx) const {
    assert(indx >= 0 && indx < bonds.size());
    return bonds.at(indx);
  }

  int resid::obBondNo(bound* bond) const {
    int post = find(bonds.begin(), bonds.end(), bond) - bonds.begin();
    if(post < bonds.size())
      return (post + 1);
    else
      return MAXNBOND;
  }

  const angle* resid::obAngl(const int indx) const {
    assert(indx >= 0 && indx < angls.size());
    return angls.at(indx);
  }

  int resid::obAnglNo(angle* angl) const {
    int post = find(angls.begin(), angls.end(), angl) - angls.begin();
    if(post < angls.size())
      return (post + 1);
    else
      return MAXNANGL;
  }

  const dihed* resid::obDhed(const int indx) const {
    assert(indx >= 0 && indx < dheds.size());
    return dheds.at(indx);
  }

  int resid::obDhedNo(dihed* dhed) const {
    int post = find(dheds.begin(), dheds.end(), dhed) - dheds.begin();
    if(post < dheds.size())
      return (post + 1);
    else
      return MAXNDHED;
  }


  myCoor resid::obCoor(const int indx) const {
    myCoor coor = 0, resM = 0;
    for(int i = 0; i < parts.size(); ++i) {
      coor += (*(parts.at(i))).obCoor(indx) * (*(parts.at(i))).obMass();
      resM += (*(parts.at(i))).obMass();
    }
    return (coor / resM);
  }

}

