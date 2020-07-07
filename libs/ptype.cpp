#include "ptype.h"
#include <cassert>

using namespace std;

namespace MD {
  //CONSTRUCTORS
  ptype::ptype() {
    vale = mass = chrg = epsl = rmin = sigm = 0;
    name = emptyNme;
  }

  ptype::ptype(const int sVale, const myCoor sMass, const myCoor sChrg,
               const myCoor sEpsl, const myCoor sRmin,
               const myCoor sSigm, const string sName) {
    vale = sVale;
    mass = sMass;
    chrg = sChrg;
    epsl = sEpsl;
    rmin = sRmin;
    sigm = sSigm;
    name = sName;
  }

  ptype::ptype(istream& iput) {
    myCoor masA;
    iput >> name >> vale >> mass >> masA >> chrg >> epsl >> rmin >> sigm;
  }

  //MODIFIERS
  const ptype& ptype::operator =(const ptype& sType) {
    if(this != &sType) {
      vale = sType.obVale();
      mass = sType.obMass();
      chrg = sType.obChrg();
      epsl = sType.obEpsl();
      rmin = sType.obRmin();
      sigm = sType.obSigm();
      name = sType.obName();
    }
    return *this;
  }
}
