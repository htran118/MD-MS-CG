/****************************************************************************
 * Angle class
 ***************************************************************************/

#ifndef ANGLE_H
#define ANGLE_H

#include "const.h"
#include "atype.h"
#include "partk.h"
#include "bound.h"
#include "resid.h"
#include "world.h"

namespace MD
{
  class atype;
  class world;
  class partk;
  class bound;
  class resid;
  class angle
  {
    public:
    //CONSTRUCTORS
      explicit angle() { size = 0; }
      angle(const partk* atom1, const partk* atom2, const partk* atom3,
            const atype* sType);
      angle(istream& iput, const world* syst);
      ~angle() { return; }

    //MODIFIERS
      //BUG: need update for bRig, bLef
      void mkLeft(const partk* sAtom) { left = sAtom; }
      void mkMidl(const partk* sAtom) { midl = sAtom; }
      void mkRigh(const partk* sAtom) { righ = sAtom; }
      void mkRsid(const resid* sMole) { mole = sMole; }
      void mkType(const atype* sType) { type = sType; }
      void mkSize(const myCoor sSize) { size = sSize; }

      const angle& operator =(const angle& angl);

    //OBSERVERS
      const partk* obLeft() const { return left; }
      const partk* obMidl() const { return midl; }
      const partk* obRigh() const { return righ; }
      const bound* obBlef() const { return bLef; }
      const bound* obBrig() const { return bRig; }
      const resid* obRsid() const { return mole; }
      const atype* obType() const { return type; }

      myCoor obRelx() const { return (*type).obRelx(); }
      myCoor obRigd() const { return (*type).obRigd(); }
      myCoor obSize() const { return size; }
      myCoor obEner() const;

      bool operator ==(const angle& angl) const;
    private:
      const partk *left, *midl, *righ;
      const bound *bLef, *bRig; 
      const resid* mole;
      const atype* type;
      myCoor size;       //relaxation length and stiffness
  };
}
#endif
