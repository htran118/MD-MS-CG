/****************************************************************************
 * Particle class
 ***************************************************************************/

#ifndef PARTK_H
#define PARTK_H

#include <vector>
#include "const.h"
#include "ptype.h"
#include "bound.h"
#include "angle.h"
#include "dihed.h"
#include "resid.h"
#include "memor.h"
#include "world.h"

namespace MD
{
  class ptype;
  class world;
  class bound;
  class angle;
  class dihed;
  class resid;
  class memor;
  class partk
  {
    public:
    //CONSTRUCTORS
      explicit partk();
      partk(const string sMark, const ptype* type);
      partk(const partk& sAtom);
      partk(istream& iput, const world* syst);
      ~partk();
    //MODIFIERS
      const partk& operator =(const partk& sAtom);

      void mkPost(const myVect* sPost) { myVectCpy(post, sPost); }
      void mkVelo(const myVect* sVelo) { myVectCpy(velo, sVelo); }
      void mkAccl(const myVect* sAccl) { myVectCpy(accl, sAccl); }
      void mkForc(const myVect* sForc);

      void mkXMem(const myVect* sPost, const int indx)
      { xMemo.mkMemo(sPost, indx); }
      void mkVMem(const myVect* sVelo, const int indx)
      { vMemo.mkMemo(sVelo, indx); }
      void mkAMem(const myVect* sAccl, const int indx)
      { aMemo.mkMemo(sAccl, indx); }
      void nwXMem(const myVect* sPost) { xMemo.nwMemo(sPost); }
      void nwVMem(const myVect* sVelo) { vMemo.nwMemo(sVelo); }
      void nwAMem(const myVect* sAccl) { aMemo.nwMemo(sAccl); }

      void mkPtyp(const ptype* sType) { type = sType; }
      void mkName(const string sName) { name = sName; }
      void mkMark(const string sMark) { mark = sMark; }
      void mkRsid(const resid* sMole) { mole = sMole; }

      void adPost(const myVect* xDif) { myVectAdd(post, xDif); }
      void adVelo(const myVect* vDif) { myVectAdd(velo, vDif); }
      void adAccl(const myVect* aDif) { myVectAdd(accl, aDif); }
      void adForc(const myVect* fDif);

      void adXMem(const int size) { xMemo.adMemo(size); }
      void adVMem(const int size) { vMemo.adMemo(size); }
      void adAMem(const int size) { aMemo.adMemo(size); }

      void adBond(bound* bond);
      void rmBond(bound* bond);
      void rmBond(const int indx);

      void adAngl(angle* angl);
      void rmAngl(angle* angl);
      void rmAngl(const int indx);

      void adDhed(dihed* dhed);
      void rmDhed(dihed* dhed);
      void rmDhed(const int indx);

      void iniTrj(istream& iput);

      void verlet();
      void lpFrgA();
      void lpFrgB();
      void vVerlA();
      void vVerlB();

      void rflect(const int indx);

    //OBSERVERS
      void wrTraj(ostream& oput) const;
      void wrCoor(ostream& oput) const;

      myCoor obCoor(const int indx) const; 

      const myVect* obPost() const { return post; }
      const myVect* obVelo() const { return velo; }
      const myVect* obAccl() const { return accl; }

      const myVect* obXMem(const int indx) const
      { return xMemo.obMemo(indx); }
      const myVect* obVMem(const int indx) const
      { return vMemo.obMemo(indx); }
      const myVect* obAMem(const int indx) const
      { return aMemo.obMemo(indx); }

      int obNbound() const { return bonds.size(); }
      const bound* obBond(const int indx) const;
      int obBondNo(const bound* bond) const;

      int obNangle() const { return angls.size(); }
      const angle* obAngl(const int indx) const;
      int obAnglNo(const angle* angl) const;

      int obNdihed() const { return dheds.size(); }
      const dihed* obDhed(const int indx) const;
      int obDhedNo(const dihed* dhed) const;

      const resid* obRsid() const { return mole; }
      const ptype* obType() const { return type; }
      string obName() const { return name; }
      string obMark() const { return mark; }
      int obVale() const { return (*type).obVale(); }
      myCoor obMass() const { return (*type).obMass(); }
      myCoor obChrg() const { return (*type).obChrg(); }
      myCoor obEpsl() const { return (*type).obEpsl(); }
      myCoor obRmin() const { return (*type).obRmin(); }
      myCoor obSigm() const { return (*type).obSigm(); }

      bool operator ==(const partk& sAtom) const
      { return (name == sAtom.obName()); }

    private:
      myVect *post, *velo, *accl, *posA, *velA, *accA;
      vector<bound*> bonds;
      vector<angle*> angls;
      vector<dihed*> dheds;
      const resid* mole;
      const ptype* type;
      string name, mark;       //name and ID number
      memor vMemo, xMemo, aMemo;
                               //memory kernel for MZ
    };

  //NON-MEMBER FUNCTIONS
  void obXdif(const partk* atom1, const partk* atom2,
              myVect* const xDif);
  void obVdif(const partk* atom1, const partk* atom2,
              myVect* const vDif);
  const bool samRsd(const partk* atom1, const partk* atom2);
  const bool samTop(const partk* atom1, const partk* atom2);
  const bool samBnd(const partk* atom1, const partk* atom2);
  const bool samAng(const partk* atom1, const partk* atom2);
  const bool samDih(const partk* atom1, const partk* atom2);
  const bound* obBond(const partk* atom1, const partk* atom2);
  const angle* obAngl(const partk* atom1, const partk* atom2,
                      const partk* atom3);
  const dihed* obDhed(const partk* atom1, const partk* atom2,
                      const partk* atom3, const partk* atom4);
}
#endif

