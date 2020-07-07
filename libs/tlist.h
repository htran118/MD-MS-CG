/****************************************************************************
 * Type-list class
 ***************************************************************************/

#ifndef UNIVERSE_H
#define UNIVERSE_H

#include <vector>
#include "rtype.h"
#include "ptype.h"
#include "btype.h"
#include "atype.h"
#include "dtype.h"

namespace MD
{
  class rtype;
  class ptype;
  class btype;
  class atype;
  class dtype;
  class tlist
  {
    public:
    //CONSTRUCTORS
      explicit tlist() { return; }
      ~tlist() { return; }

    //MODIFIERS
      void adType(istream& iput);
      void adRtyp(istream& iput);
      void adPtyp(istream& iput);
      void adBtyp(istream& iput);
      void adAtyp(istream& iput);
      void adDtyp(istream& iput);

    //OBSERVERS
      const vector<rtype>& obRtypLs() const { return rtyps; }
      const vector<ptype>& obPtypLs() const { return ptyps; }
      const vector<btype>& obBtypLs() const { return btyps; }
      const vector<atype>& obAtypLs() const { return atyps; }
      const vector<dtype>& obDtypLs() const { return dtyps; }

      const rtype* obRtyp(const string name) const;
      const ptype* obPtyp(const string name) const;
      const btype* obBtyp(const ptype* type1, const ptype* type2) const;
      const atype* obAtyp(const ptype* type1, const ptype* type2,
                          const ptype* type3) const;
      const dtype* obDtyp(const ptype* type1, const ptype* type2,
                          const ptype* type3, const ptype* type4) const;

      int obRtypNo(const rtype* rtyp) const;
      int obPtypNo(const ptype* ptyp) const;
      int obBtypNo(const btype* btyp) const;
      int obAtypNo(const atype* atyp) const;
      int obDtypNo(const dtype* dtyp) const;

      const rtype* obRtyp(const int indx) const;
      const ptype* obPtyp(const int indx) const;
      const btype* obBtyp(const int indx) const;
      const atype* obAtyp(const int indx) const;
      const dtype* obDtyp(const int indx) const;

      int obNrtype() const { return rtyps.size(); }
      int obNptype() const { return ptyps.size(); }
      int obNbtype() const { return btyps.size(); }
      int obNatype() const { return atyps.size(); }
      int obNdtype() const { return dtyps.size(); }

    private:
      vector<rtype> rtyps;
      vector<ptype> ptyps;
      vector<btype> btyps;
      vector<atype> atyps;
      vector<dtype> dtyps;
  };
}
#endif
