/****************************************************************************
 * Residue class
 ***************************************************************************/

#ifndef RESID_H
#define RESID_H

#include <vector>
#include "partk.h"
#include "rtype.h"
#include "bound.h"
#include "angle.h"
#include "dihed.h"
#include "tlist.h"

namespace MD
{
  class rtype;
  class tlist;
  class partk;
  class bound;
  class angle;
  class dihed;
  class resid
  {
    public:
    //CONSTRUCTORS
      explicit resid() { return; }
      resid(istream& iput, const tlist* types);
      ~resid() { return; }
    //MODIFIERS
      void adPart(partk* atom);
      void adBond(bound* bond);
      void adAngl(angle* angl);
      void adDhed(dihed* dhed);

      void rmPart(partk* atom);
      void rmPart(const int indx);
      void rmBond(bound* bond);
      void rmBond(const int indx);
      void rmAngl(angle* angl);
      void rmAngl(const int indx);
      void rmDhed(dihed* dhed);
      void rmDhed(const int indx);

      const resid& operator =(const resid& mole);

      void rflect(const int indx);

    //OBSERVERS
      const partk* obPart(const int indx) const;
      const bound* obBond(const int indx) const;
      const angle* obAngl(const int indx) const;
      const dihed* obDhed(const int indx) const;
      int obPartNo(partk* atom) const;
      int obBondNo(bound* bond) const;
      int obAnglNo(angle* angl) const;
      int obDhedNo(dihed* dhed) const;

      int obNpartk() const { return parts.size(); }
      int obNbound() const { return bonds.size(); }
      int obNangle() const { return angls.size(); }
      int obNdihed() const { return dheds.size(); }

      myCoor obCoor(const int indx) const;

      string obName() const { return name; }
      const rtype* obRtyp() const { return type; }

    private:
      vector<partk*> parts;
      vector<bound*> bonds;
      vector<angle*> angls;
      vector<dihed*> dheds;
      string name;
      const rtype *type;
  };
}
#endif
