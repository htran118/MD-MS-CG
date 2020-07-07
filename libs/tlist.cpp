#include <cassert>
#include <sstream>
#include <fstream>
#include <algorithm>
#include "tlist.h"

using namespace std;

namespace MD {
  //CONSTRUCTORS

  //MODIFIERS
  void tlist::adType(istream& iput) {
    string ctrl;
    while(iput >> ctrl) {
      if(ctrl == ctrlIgnr)
        iput.ignore(256, '\n');
      else if(ctrl == ctrlAtom) {
        ptyps.push_back(ptype(iput)); }
      else if(ctrl == ctrlBond)
        btyps.push_back(btype(iput, this));
      else if(ctrl == ctrlAngl)
        atyps.push_back(atype(iput, this));
      else if(ctrl == ctrlDhed)
        dtyps.push_back(dtype(iput, this));
      else if(ctrl == ctrlRsid)
        rtyps.push_back(rtype(iput, this));
    }
  }

  void tlist::adRtyp(istream& iput) {
    rtyps.push_back(rtype(iput, this));
  }

  void tlist::adPtyp(istream& iput) {
    ptyps.push_back(ptype(iput));
  }

  void tlist::adBtyp(istream& iput) {
    btyps.push_back(btype(iput, this));
  }

  void tlist::adAtyp(istream& iput) {
    atyps.push_back(atype(iput, this));
  }

  void tlist::adDtyp(istream& iput) {
    dtyps.push_back(dtype(iput, this));
  }

  //OBSERVERS
  const rtype* tlist::obRtyp(const string name) const {
    for(int i = 0; i < rtyps.size(); ++i)
      if((rtyps.at(i)).obName() == name)
        return &rtyps.at(i);
    cerr << "Cannot find resid type " << name << endl;
    return NULL;
  }

  const ptype* tlist::obPtyp(const string name) const {
    for(int i = 0; i < ptyps.size(); ++i)
      if((ptyps.at(i)).obName() == name)
        return &ptyps.at(i);
    cerr << "Cannot find partk type " << name << endl;
    return NULL;
  }

  const btype* tlist::obBtyp(const ptype* type1,
                             const ptype* type2) const {
    for(int i = 0; i < btyps.size(); ++i)
      if(((btyps.at(i)).obLeft() == type1 &&
          (btyps.at(i)).obRigh() == type2) ||
         ((btyps.at(i)).obLeft() == type2 &&
          (btyps.at(i)).obRigh() == type1))
        return &btyps.at(i);
    cerr << "Cannot find bond type connecting " << (*type1).obName()
         << " to " << (*type2).obName() << endl;
    return NULL;
  }

  const atype* tlist::obAtyp(const ptype* type1, const ptype* type2,
                             const ptype* type3) const {
    for(int i = 0; i < atyps.size(); ++i)
      if(((atyps.at(i)).obLeft() == type1 &&
          (atyps.at(i)).obMidl() == type2 &&
          (atyps.at(i)).obRigh() == type3) ||
         ((atyps.at(i)).obLeft() == type3 &&
          (atyps.at(i)).obMidl() == type2 &&
          (atyps.at(i)).obRigh() == type1))
        return &atyps.at(i);
    cerr << "Cannot find angle type connecting " << (*type2).obName()
         << " to " << (*type1).obName()
         << " and " << (*type3).obName() << endl;
    return NULL;
  }

  const dtype* tlist::obDtyp(const ptype* type1, const ptype* type2,
                             const ptype* type3, const ptype* type4) const {
    for(int i = 0; i < dtyps.size(); ++i)
      if(((dtyps.at(i)).obUppr() == type1 &&
          (dtyps.at(i)).obLeft() == type2 &&
          (dtyps.at(i)).obRigh() == type3 &&
          (dtyps.at(i)).obLowr() == type4) ||
         ((dtyps.at(i)).obUppr() == type4 &&
          (dtyps.at(i)).obLeft() == type3 &&
          (dtyps.at(i)).obRigh() == type2 &&
          (dtyps.at(i)).obLowr() == type1))
        return &dtyps.at(i);
    cerr << "Cannot find dihedral type connecting " << (*type2).obName()
         << " and " << (*type3).obName()
         << " to " << (*type1).obName()
         << " and " << (*type4).obName() << endl;
    return NULL;
  }

  const rtype* tlist::obRtyp(const int indx) const {
    assert(indx >= 0 && indx < rtyps.size());
    return &rtyps.at(indx);
  }

  const ptype* tlist::obPtyp(const int indx) const {
    assert(indx >= 0 && indx < ptyps.size());
    return &ptyps.at(indx);
  }

  const btype* tlist::obBtyp(const int indx) const {
    assert(indx >= 0 && indx < btyps.size());
    return &btyps.at(indx);
  }

  const atype* tlist::obAtyp(const int indx) const {
    assert(indx >= 0 && indx < atyps.size());
    return &atyps.at(indx);
  }

  const dtype* tlist::obDtyp(const int indx) const {
    assert(indx >= 0 && indx < dtyps.size());
    return &dtyps.at(indx);
  }
}
