#include "btype.h"
#include <cassert>

using namespace std;

namespace MD {
  //CONSTRUCTORS
  btype::btype(const ptype *type1, const ptype* type2, const myCoor sRelx,
               const myCoor sRigd) {
    relx = sRelx;
    rigd = sRigd;
    left = type1;
    righ = type2;
  }

  btype::btype(istream& iput, const tlist* types) {
    string name1, name2;
    iput >> name1 >> name2 >> rigd >> relx;
    left =(*types).obPtyp(name1);
    righ =(*types).obPtyp(name2);
  }
      
  //MODIFIERS
  const btype& btype::operator =(const btype& sType) {
    if(this != &sType) {
      relx = sType.obRelx();
      rigd = sType.obRigd();
      left = sType.obLeft();
      righ = sType.obRigh();
    }
    return *this;
  }

  bool btype::operator ==(const btype& sType) const {
    return (left == sType.obLeft() && righ == sType.obRigh());
  }
}
