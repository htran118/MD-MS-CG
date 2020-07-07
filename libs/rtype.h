/****************************************************************************
 * Residue-type class
 ***************************************************************************/

#ifndef RTYPE_H
#define RTYPE_H
#include <map>
#include <vector>
#include "tlist.h"
#include "ptype.h"
#include "btype.h"
#include "atype.h"
#include "dtype.h"

using namespace std;

namespace MD
{
  class tlist;
  class ptype;
  class btype;
  class atype;
  class dtype;
  class doubl
  {
    public:
    //CONSTRUCTORS
      explicit doubl();
      doubl(const string name1, const string name2);
      doubl(const doubl& group);
      doubl(istream& iput);
    //MODIFIERS
      const doubl& operator =(const doubl& group);
    //OBSERVERS
      string obLeft() const { return leftNm; }
      string obRigh() const { return righNm; }
      bool operator ==(const doubl& group) const
      { return (leftNm == group.obLeft() && righNm == group.obRigh()); }
      bool operator <(const doubl& group) const;
    private:
      string leftNm, righNm;
  };

  class tripl
  {
    public:
    //CONSTRUCTORS
      explicit tripl();
      tripl(const string name1, const string name2, const string name3);
      tripl(const tripl& group);
      tripl(istream& iput);
    //MODIFIERS
      const tripl& operator =(const tripl& group);
    //OBSERVERS
      string obLeft() const { return leftNm; }
      string obMidl() const { return midlNm; }
      string obRigh() const { return righNm; }
      bool operator ==(const tripl& group) const
      { return (leftNm == group.obLeft() && midlNm == group.obMidl() &&
                righNm == group.obRigh()); }
      bool operator <(const tripl& group) const;
    private:
      string leftNm, midlNm, righNm;
  };

  class quadr
  {
     public:
    //CONSTRUCTORS
      explicit quadr();
      quadr(const string name1, const string name2,
            const string name3, const string name4);
      quadr(const quadr& group);
      quadr(istream& iput);
    //MODIFIERS
      const quadr& operator =(const quadr& group);
    //OBSERVERS
      string obUppr() const { return upprNm; }
      string obLeft() const { return leftNm; }
      string obRigh() const { return righNm; }
      string obLowr() const { return lowrNm; }
      bool operator ==(const quadr& group) const
      { return (upprNm == group.obUppr() && leftNm == group.obLeft() &&
                righNm == group.obRigh() && lowrNm == group.obLowr()); }
      bool operator <(const quadr& group) const;
    private:
      string upprNm, leftNm, righNm, lowrNm;
  };

  class rtype
  {
    public:
    //CONSTRUCTORS
      explicit rtype();
      rtype(istream& iput, const tlist* types);
      ~rtype();
    //MODIFIERS
      void adPtypMap(istream& iput, const tlist* types);
      void adBtypMap(istream& iput, const tlist* types);
      void adAtypMap(istream& iput, const tlist* types);
      void adDtypMap(istream& iput, const tlist* types);
      const rtype& operator =(const rtype& sType);
    //OBSERVERS
      string obName() const { return name; }

      const map <string, const ptype*>& obPtypMap() const
      { return partToPtyp; }
      const map <doubl, const btype*>& obBtypMap() const
      { return partToBtyp; }
      const map <tripl, const atype*>& obAtypMap() const
      { return partToAtyp; }
      const map <quadr, const dtype*>& obDtypMap() const
      { return partToDtyp; }

      const ptype* obPtyp(const string sName) const
      { return partToPtyp.at(sName); }
      const btype* obBtyp(const string name1, const string name2) const 
      { return partToBtyp.at(doubl(name1, name2)); }
      const atype* obAtyp(const string name1, const string name2,
                          const string name3) const
      { return partToAtyp.at(tripl(name1, name2, name3)); }
      const dtype* obDtyp(const string name1, const string name2,
                          const string name3, const string name4) const
      { return partToDtyp.at(quadr(name1, name2, name3, name4)); }

      int obNptype() const { return partToPtyp.size(); }

      bool operator ==(const rtype& sType) const
      { return (name == sType.obName()); }
    private:
      map <string, const ptype*> partToPtyp;
      map <doubl, const btype*> partToBtyp;
      map <tripl, const atype*> partToAtyp;
      map <quadr, const dtype*> partToDtyp;
      string name;
  };
}
#endif
