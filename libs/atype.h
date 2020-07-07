/****************************************************************************
 * Angle-type class
 ***************************************************************************/

#ifndef ATYPE_H
#define ATYPE_H

#include <vector>
#include "ptype.h"
#include "tlist.h"

namespace MD
{
  class tlist;
  class ptype;
  class atype
  {
    public:
    //CONSTRUCTORS
      explicit atype() { relx = rigd = 0.0; }
      atype(const ptype* type1, const ptype* type2, const ptype* type3,
            const myCoor sRelx, const myCoor sRigd);
      atype(istream& iput, const tlist* types);
      ~atype() { return; }
    //MODIFIERS
      const atype& operator =(const atype& sType);
    //OBSERVERS
      myCoor obRelx() const { return relx; }
      myCoor obRigd() const { return rigd; }
      const ptype* obMidl() const { return midl; }
      const ptype* obLeft() const { return left; }
      const ptype* obRigh() const { return righ; }
      bool operator ==(const atype& sType) const;
    private:
      const ptype *midl, *left, *righ;
      myCoor relx, rigd;
  };
}
#endif
