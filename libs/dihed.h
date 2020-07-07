/****************************************************************************
 * Dihed class
 ***************************************************************************/

#ifndef DIHED_H
#define DIHED_H

#include "const.h"
#include "dtype.h"
#include "partk.h"
#include "bound.h"
#include "angle.h"
#include "resid.h"
#include "world.h"

namespace MD
{
  class dtype;
  class world;
  class partk;
  class bound;
  class angle;
  class resid;
  class dihed
  {
    public:
    //CONSTRUCTORS
      explicit dihed() { size = 0; }
      dihed(const partk* atom1, const partk* atom2, const partk* atom3,
            const partk* atom4, const dtype* sType);
      dihed(istream& iput, const world* syst);
      ~dihed() { return; }

    //MODIFIERS
      //BUG: need update for bRig, bMid, bLef, aLef, aRig
      void mkUppr(const partk* sAtom) { uppr = sAtom; }
      void mkLeft(const partk* sAtom) { left = sAtom; }
      void mkRigh(const partk* sAtom) { righ = sAtom; }
      void mkLowr(const partk* sAtom) { lowr = sAtom; }
      void mkRsid(const resid* sMole) { mole = sMole; }
      void mkType(const dtype* sType) { type = sType; }
      void mkSize(const myCoor sSize) { size = sSize; }

      const dihed& operator =(const dihed& dhed);

    //OBSERVERS
      const partk* obUppr() const { return uppr; }
      const partk* obLeft() const { return left; }
      const partk* obRigh() const { return righ; }
      const partk* obLowr() const { return lowr; }
      const bound* obBlef() const { return bLef; }
      const bound* obBmid() const { return bMid; }
      const bound* obBrig() const { return bRig; }
      const angle* obAlef() const { return aLef; }
      const angle* obArig() const { return aRig; }
      const resid* obRsid() const { return mole; }
      const dtype* obType() const { return type; }

      myCoor obRelx() const { return (*type).obRelx(); }
      myCoor obRig1() const { return (*type).obRig1(); }
      myCoor obRig2() const { return (*type).obRig2(); }
      myCoor obRig3() const { return (*type).obRig3(); }
      myCoor obSize() const { return size; }
      myCoor obEner() const;

      bool operator ==(const dihed& dhed) const;
    private:
      const partk *uppr, *left, *righ, *lowr;
      const bound *bLef, *bMid, *bRig;
      const angle *aLef, *aRig;
      const resid* mole;
      const dtype* type;
      myCoor size;                //relaxation length
  };
}
#endif
