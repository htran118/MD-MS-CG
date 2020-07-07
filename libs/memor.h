/********************************************************************
 * Container class for the memory
 *******************************************************************/

#ifndef MEMOR_H
#define MEMOR_H

#include <vector>
#include "const.h"

namespace MD
{
  class partk;
  class memor
  {
    public:
    //CONSTRUCTORS
      explicit memor() { mPos = 0; }
      ~memor();
    //MODIFIERS
      const memor& operator =(const memor& sMemo);
      void adMemo(const int size);
      void fitSze() { memos.shrink_to_fit(); }
      void nwMemo(const myVect* sVect);
      void mkMemo(const myVect* sVect, const int indx);
    //OBSERVERS
      const myVect* obMemo(const int indx) const;
      int obNmemor() const { return memos.size(); }
      int obMpos() const { return mPos; }
    protected:
      vector<myVect*> memos;
      int mPos;
  };

  //container class for pair memory (just a fancy 2D->1D array)
  //data for same pair at different time are contigious
  class pmemo: public memor
  {
    public:
    //CONSTRUCTORS
      explicit pmemo();
      ~pmemo() { return; }
    //MODIFIERS
      template <class type>
      void adPmem(const int size, const vector<type>* objcs);
      //void adPmem(const int size, const vector<partk>* objcs);
      //NOTE: nwPmem does not increment the current pointer, so we need 
      //nwMpos to update
      template <class type>
      void nwPmem(const myVect* sVect, const vector<type>* objcs,
                  const type* obj1, const type* obj2);
      //void nwPmem(const myVect* sVect, const vector<partk>* objcs,
      //            const partk* obj1, const partk* obj2);
      void nwMpos();
      template <class type>
      void mkPmem(const myVect* sVect, const vector<type>* parts,
                  const type* atom1, const type* atom2, const int indx);
      //void mkPmem(const myVect* sVect, const vector<partk>* parts,
      //            const partk* atom1, const partk* atom2, const int indx);
    //OBSERVERS
      template <class type>
      const myVect* obPmem(const vector<type>* parts, const type* atom1,
                           const type* atom2, const int indx) const;
      //const myVect* obPmem(const vector<partk>* parts, const partk* atom1,
      //                     const partk* atom2, const int indx) const;
      int obNpmemo() const { return mLen; }
    private:
      int mLen;
      bool updt;
  };

  //container for auxiliary variables
  class auxMm
  {
    public:
    //CONSTRUCTORS
      explicit auxMm() { return; }
      auxMm(const auxMm &sMemo);
      auxMm(const int size);
      ~auxMm();
    //MODIFIERS
      const auxMm& operator =(const auxMm& sMemo);
      void adMemo(const int size);
      void iniTrj(istream& iput);
      void mkPost(const myVect* sPost) { myVectCpy(post, sPost); }
      void mkVelo(const myVect* sVelo) { myVectCpy(velo, sVelo); }
      void adPost(const myVect* xDif) { myVectAdd(post, xDif); }
      void adVelo(const myVect* vDif) { myVectAdd(velo, vDif); }
      void nxtStp();
    //OBSERVERS
      void wrTraj(ostream& oput) const;
      const myVect* obPost() const { return post; }
      const myVect* obVelo() const { return velo; }
    private:
      myVect *post, *velo;
  };
}

#include "partk.h"

#endif
