/****************************************************************************
 * Bond-type class
 ***************************************************************************/

#ifndef BTYPE_H
#define BTYPE_H

#include <vector>
#include "ptype.h"
#include "tlist.h"

namespace MD
{
  class tlist;
  class ptype;
  class btype
  {
    public:
    //CONSTRUCTORS
      explicit btype() { relx = rigd = 0.0; }
      btype(const ptype* type1, const ptype* type2, const myCoor sRelx,
            const myCoor sRigd);
      btype(istream& iput, const tlist* types);
      ~btype() { return; }
    //MODIFIERS
      const btype& operator =(const btype& sType);
    //OBSERVERS
      myCoor obRelx() const { return relx; }
      myCoor obRigd() const { return rigd; }
      const ptype* obLeft() const { return left; }
      const ptype* obRigh() const { return righ; }
      bool operator ==(const btype& sType) const;
    private:
      const ptype *left, *righ;
      myCoor relx, rigd;
  };
}
#endif
