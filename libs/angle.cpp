#include "angle.h"
#include <cassert>
#include <cmath>

using namespace std;

namespace MD
{
  //CONSTRUCTORS

  angle::angle(const partk* atom1, const partk* atom2, const partk* atom3,
               const atype* sType) {
    assert(!(atom1 == atom2 || atom2 == atom3 || atom3 == atom1) &&
            (*atom1).obRsid() == (*atom2).obRsid() &&
            (*atom3).obRsid() == (*atom2).obRsid());
    left = atom1;
    midl = atom2;
    righ = atom3;
    bLef = obBond(left, midl);
    bRig = obBond(righ, midl);
    mole = (*midl).obRsid();
    type = sType;
  }

  angle::angle(istream& iput, const world* syst) {
    string name1, name2, name3;
    iput >> name1 >> name2 >> name3;
    left = (*syst).obPart(name1);
    midl = (*syst).obPart(name2);
    righ = (*syst).obPart(name3);
    if(left != NULL && midl != NULL && righ != NULL) {
      assert((*left).obRsid() == (*midl).obRsid() &&
             (*righ).obRsid() == (*midl).obRsid());
      mole = (*midl).obRsid();
      bLef = obBond(left, midl);
      bRig = obBond(righ, midl);
      type = (*((*mole).obRtyp())).obAtyp((*left).obMark(), (*midl).obMark(),
                                          (*righ).obMark());
    }
  }

  //MODIFIERS
  const angle& angle::operator =(const angle& angl) {
    if(this != &angl) {
      type = angl.obType();
    }
    return *this;
  }

  //OBSERVERS
  myCoor angle::obEner() const { 
    return (0.5 * (*type).obRigd() * (size - (*type).obRelx())
                                   * (size - (*type).obRelx()));
  }

  bool angle::operator ==(const angle& angl) const {
    return (left == angl.obLeft() && midl == angl.obMidl() &&
            righ == angl.obRigh());
  }
}
