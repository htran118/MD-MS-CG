#include "atype.h"
#include <cassert>

using namespace std;

namespace MD {
  //CONSTRUCTORS
  atype::atype(const ptype* type1, const ptype* type2, const ptype* type3,
               const myCoor sRelx, const myCoor sRigd) {
    relx = sRelx;
    rigd = sRigd;
    left = type1;
    midl = type2;
    righ = type3;
  }

  atype::atype(istream& iput, const tlist* types) {
    string name1, name2, name3;
    iput >> name1 >> name2 >> name3 >> rigd >> relx;
    relx *= ONEPIE / 180.0;
    left =(*types).obPtyp(name1);
    midl =(*types).obPtyp(name2);
    righ =(*types).obPtyp(name3);
  }

  //MODIFIERS
  const atype& atype::operator =(const atype& sType) {
    if(this != &sType) {
      relx = sType.obRelx();
      rigd = sType.obRigd();
      left = sType.obLeft();
      midl = sType.obMidl();
      righ = sType.obRigh();
    }
    return *this;
  }

  bool atype::operator ==(const atype& sType) const {
    return (righ == sType.obRigh() && left == sType.obLeft() &&
            midl == sType.obMidl());
  }
}
