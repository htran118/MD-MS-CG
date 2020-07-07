/****************************************************************************
 * Bond class
 ***************************************************************************/

#ifndef HBOND_H
#define HBOND_H

#include <vector>
#include "const.h"
#include "btype.h"
#include "partk.h"
#include "resid.h"
#include "world.h"

namespace MD
{
  class btype;
  class world;
  class partk;
  class resid;
  class bound
  {
    public:
    //CONSTRUCTORS
      explicit bound() { size = 0; }
      bound(const partk* atom1, const partk* atom2, const btype* sType);
      bound(istream& iput, const world* syst);
      ~bound() { return; }

    //MODIFIERS
      void mkLeft(const partk* sAtom) { left = sAtom; }
      void mkRigh(const partk* sAtom) { righ = sAtom; }
      void mkRsid(const resid* sMole) { mole = sMole; }
      void mkType(const btype* sType) { type = sType; }
      void mkSize(const myCoor sSize) { size = sSize; }

      const bound& operator =(const bound& bond);

    //OBSERVERS
      const partk* obLeft() const { return left; }
      const partk* obRigh() const { return righ; }
      const resid* obRsid() const { return mole; }
      const btype* obType() const { return type; }

      myCoor obRelx() const { return (*type).obRelx(); }
      myCoor obRigd() const { return (*type).obRigd(); }
      myCoor obSize() const { return size; }
      myCoor obEner() const;

      bool operator ==(const bound& bond) const;
    private:
      const partk *left, *righ;
      const resid* mole;
      const btype* type;
      myCoor size;         //relaxation length and stiffness
  };
}
#endif
