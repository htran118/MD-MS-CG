#include "rtype.h"

using namespace std;

namespace MD {
  //class DOUBLET
  doubl::doubl() {
    leftNm = righNm = emptyNme;
  }

  doubl::doubl(const string name1, const string name2) {
    leftNm = name1;
    righNm = name2;
  }

  doubl::doubl(const doubl& group) {
    leftNm = group.obLeft();
    righNm = group.obRigh();
  }

  doubl::doubl(istream& iput) {
    iput >> leftNm >> righNm;
  }

  const doubl& doubl::operator =(const doubl& group) {
    leftNm = group.obLeft();
    righNm = group.obRigh();
    return *this;
  }

  bool doubl::operator <(const doubl& group) const {
    if(leftNm == group.obLeft())
      return (righNm < group.obRigh());
    else
      return (leftNm < group.obLeft());
  }

  //class TRIPLET
  tripl::tripl() {
    leftNm = midlNm = righNm = emptyNme;
  }

  tripl::tripl(const string name1, const string name2, const string name3) {
    leftNm = name1;
    midlNm = name2;
    righNm = name3;
  }

  tripl::tripl(const tripl& group) {
    leftNm = group.obLeft();
    midlNm = group.obMidl();
    righNm = group.obRigh();
  }

  tripl::tripl(istream& iput) {
    iput >> leftNm >> midlNm >> righNm;
  }

  const tripl& tripl::operator =(const tripl& group) {
    leftNm = group.obLeft();
    midlNm = group.obMidl();
    righNm = group.obRigh();
    return *this;
  }

  bool tripl::operator <(const tripl& group) const {
    if(leftNm == group.obLeft())
      if(midlNm == group.obMidl())
        return (righNm < group.obRigh());
      else
        return (midlNm < group.obMidl());
    else
      return (leftNm < group.obLeft());
  }

  //class QUADRUPLET
  quadr::quadr() {
    upprNm = leftNm = righNm = lowrNm =  emptyNme;
  }

  quadr::quadr(const string name1, const string name2,
               const string name3, const string name4) {
    upprNm = name1;
    leftNm = name2;
    righNm = name3;
    lowrNm = name4;
  }

  quadr::quadr(const quadr& group) {
    upprNm = group.obUppr();
    leftNm = group.obLeft();
    righNm = group.obRigh();
    lowrNm = group.obLowr();
  }

  quadr::quadr(istream& iput) {
    iput >> upprNm >> leftNm >> righNm >> lowrNm;
  }

  const quadr& quadr::operator =(const quadr& group) {
    upprNm = group.obUppr();
    leftNm = group.obLeft();
    righNm = group.obRigh();
    lowrNm = group.obLowr();
    return *this;
  }

  bool quadr::operator <(const quadr& group) const {
    if(upprNm == group.obUppr())
      if(leftNm == group.obLeft())
        if(righNm == group.obRigh())
          return (lowrNm < group.obLowr());
        else
          return (righNm < group.obRigh());
      else
        return (leftNm < group.obLeft());
    else
      return (upprNm < group.obUppr());
  }
 
  //class RTYPE
  rtype::rtype() {
    name = emptyNme;
  }

  rtype::~rtype() {
    partToPtyp.clear();
    partToBtyp.clear();
    partToAtyp.clear();
    partToDtyp.clear();
  }

  rtype::rtype(istream& iput, const tlist* types) {
    string ctrl;
    iput >> name;
    while(ctrl != ctrlStop) {
      iput >> ctrl;
      if(ctrl == ctrlIgnr)
        iput.ignore(256, '\n');
      else if(ctrl == ctrlAtom)
        this->adPtypMap(iput, types);
      else if(ctrl == ctrlBond)
        this->adBtypMap(iput, types);
      else if(ctrl == ctrlAngl)
        this->adAtypMap(iput, types);
      else if(ctrl == ctrlDhed)
        this->adDtypMap(iput, types);
    }
  }

  const rtype& rtype::operator =(const rtype& sType) {
    if(this != &sType) {
      name = sType.obName();
      partToPtyp = sType.obPtypMap();
      partToBtyp = sType.obBtypMap();
      partToAtyp = sType.obAtypMap();
      partToDtyp = sType.obDtypMap();
    }
    return *this;
  }

  inline void rtype::adPtypMap(istream& iput, const tlist* types) {
    string partNm, ptypeNm;
    iput >> partNm >> ptypeNm;
    const ptype* sType = (*types).obPtyp(ptypeNm);
    if(sType != NULL)
      partToPtyp.insert(pair <string, const ptype*>(partNm, sType));
  }

  inline void rtype::adBtypMap(istream& iput, const tlist* types) {
    doubl group = doubl(iput);
    string name1 = group.obLeft(), name2 = group.obRigh();
    const ptype *type1 = partToPtyp.at(name1),
                *type2 = partToPtyp.at(name2);
    const btype* sType = (*types).obBtyp(type1, type2);
    if(sType != NULL)
      partToBtyp.insert(pair <doubl, const btype*>(group, sType));
  }

  inline void rtype::adAtypMap(istream& iput, const tlist* types) {
    tripl group = tripl(iput);
    string name1 = group.obLeft(), name2 = group.obMidl(),
           name3 = group.obRigh();
    const ptype *type1 = partToPtyp.at(name1),
                *type2 = partToPtyp.at(name2),
                *type3 = partToPtyp.at(name3);
    const atype* sType = (*types).obAtyp(type1, type2, type3);
    if(sType != NULL)
      partToAtyp.insert(pair <tripl, const atype*>(group, sType));
  }

  inline void rtype::adDtypMap(istream& iput, const tlist* types) {
    quadr group = quadr(iput);
    string name1 = group.obUppr(), name2 = group.obLeft(),
           name3 = group.obRigh(), name4 = group.obLowr();
    const ptype *type1 = partToPtyp.at(name1),
                *type2 = partToPtyp.at(name2),
                *type3 = partToPtyp.at(name3),
                *type4 = partToPtyp.at(name4);
    const dtype* sType = (*types).obDtyp(type1, type2, type3, type4);
    if(sType != NULL)
      partToDtyp.insert(pair <quadr, const dtype*>(group, sType));
  }
}
