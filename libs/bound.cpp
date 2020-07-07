#include "bound.h"
#include <cassert>
#include <cmath>

using namespace std;

namespace MD {
  //CONSTRUCTORS

  bound::bound(const partk* atom1, const partk* atom2, const btype* sType) {
    assert(!(atom1 == atom2) &&
           (*atom1).obRsid() == (*atom2).obRsid());
    left = atom1;
    righ = atom2;
    mole = (*left).obRsid();
    type = sType;
  }

  bound::bound(istream& iput, const world* systm) {
    string name1, name2;
    iput >> name1 >> name2;
    left = (*systm).obPart(name1);
    righ = (*systm).obPart(name2);
    if(left != NULL && righ != NULL) {
      assert((*left).obRsid() == (*righ).obRsid());
      mole = (*left).obRsid();
      type = (*((*mole).obRtyp())).obBtyp((*left).obMark(), (*righ).obMark());
    }
  }
    
  //MODIFIERS
  const bound& bound::operator =(const bound& bond) {
    if(this != &bond) {
      type = bond.obType();
    }
    return *this;
  }

  //OBSERVERS
  myCoor bound::obEner() const {
    return (0.5 * (*type).obRigd() * (size - (*type).obRelx())
                                   * (size - (*type).obRelx()));
  }

  bool bound::operator ==(const bound& bond) const {
    return (left == bond.obLeft() && righ == bond.obRigh());
  }
}
