/****************************************************************************
 * Particle-type class
 ***************************************************************************/

#ifndef PTYPE_H
#define PTYPE_H

#include "const.h"

using namespace std;

namespace MD
{
  class ptype
  {
    public:
    //CONSTRUCTORS
      explicit ptype();
      ptype(const int sVale, const myCoor sMass, const myCoor sChrg,
            const myCoor sEpsl, const myCoor sRmin,
            const myCoor sSigm, const string sName);
      ptype(istream& iput);
      ~ptype() { return; }
    //MODIFIERS
      const ptype& operator =(const ptype& sType);
    //OBSERVERS
      int obVale() const { return vale; }
      myCoor obMass() const { return mass; }
      myCoor obChrg() const { return chrg; }
      myCoor obEpsl() const { return epsl; }
      myCoor obRmin() const { return rmin; }
      myCoor obSigm() const { return sigm; }
      string obName() const { return name; }
      //for testing purpose only
      bool operator ==(const ptype& sType) const
      { return (name == sType.obName()); }
    private:
      int vale;
      myCoor mass, chrg, epsl, rmin, sigm;
      string name;
   };
}
#endif
