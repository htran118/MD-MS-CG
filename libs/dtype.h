/****************************************************************************
 * Dihedral-type class
 ***************************************************************************/

#ifndef DTYPE_H
#define DTYPE_H

#include <vector>
#include "ptype.h"
#include "tlist.h"

namespace MD
{
  class tlist;
  class ptype;
  class dtype
  {
    public:
    //CONSTRUCTORS
      explicit dtype() { relx = rig1 = rig2 = rig3 = 0.0; };
      dtype(const ptype* type1, const ptype* type2, const ptype* type3,
            const ptype* type4, const myCoor sRelx, const myCoor sRig1,
            const myCoor sRig2, const myCoor sRig3);
      dtype(istream& iput, const tlist* types);
      ~dtype() { return; }
    //MODIFIERS
      const dtype& operator =(const dtype& sType);
    //OBSERVERS
      myCoor obRelx() const { return relx; }
      myCoor obRig1() const { return rig1; }
      myCoor obRig2() const { return rig2; }
      myCoor obRig3() const { return rig3; }
      const ptype* obUppr() const { return uppr; }
      const ptype* obLeft() const { return left; }
      const ptype* obRigh() const { return righ; }
      const ptype* obLowr() const { return lowr; }
      bool operator ==(const dtype& sType) const;
    private:
      const ptype *uppr, *left, *righ, *lowr;
      myCoor relx, rig1, rig2, rig3;
  };
}
#endif
