#include "dtype.h"
#include <cassert>

using namespace std;

namespace MD {
  //CONSTRUCTORS
  dtype::dtype(const ptype* type1, const ptype* type2, const ptype* type3,
               const ptype* type4, const myCoor sRelx, const myCoor sRig1,
               const myCoor sRig2, const myCoor sRig3) {
    relx = sRelx;
    rig1 = sRig1;
    rig2 = sRig2;
    rig3 = sRig3;
    uppr = type1;
    left = type2;
    righ = type3;
    lowr = type4;
  }

  dtype::dtype(istream& iput, const tlist* types) {
    string name1, name2, name3, name4;
    iput >> name1 >> name2 >> name3 >> name4 >> rig1 >> rig2 >> rig3 >> relx;
    relx *= ONEPIE / 180.0;
    uppr =(*types).obPtyp(name1);
    left =(*types).obPtyp(name2);
    righ =(*types).obPtyp(name3);
    lowr =(*types).obPtyp(name4);
  }

  //MODIFIERS
  const dtype& dtype::operator =(const dtype& sType) {
    if(this != &sType) {
      relx = sType.obRelx();
      rig1 = sType.obRig1();
      rig2 = sType.obRig2();
      rig3 = sType.obRig3();
      uppr = sType.obUppr();
      left = sType.obLeft();
      righ = sType.obRigh();
      lowr = sType.obLowr();
    }
    return *this;
  }

  bool dtype::operator ==(const dtype& sType) const {
    return (uppr == sType.obUppr() && righ == sType.obRigh() &&
            left == sType.obLeft() && lowr == sType.obLowr());
  }
}
