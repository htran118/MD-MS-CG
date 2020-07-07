#include <cassert>
#include <climits>
#include "memor.h"

using namespace std;

namespace MD {
  //class MEMOR
  memor::~memor() {
    for(int i = 0; i < memos.size(); ++i)
      myVectMemdel(memos.at(i));
  }

  const memor& memor::operator =(const memor& sMemo) {
    if(memos.size() != sMemo.obNmemor()) {
      if(memos.size() < sMemo.obNmemor())
        for(int i = 0; i < (sMemo.obNmemor() - memos.size()); ++i)
          memos.push_back(myVectMalloc(DIM));
      if(memos.size() > sMemo.obNmemor())
        for(int i = 0; i < (memos.size() - sMemo.obNmemor()); ++i)
          memos.pop_back();
    }
    for(int i = 0; i < sMemo.obNmemor(); ++i)
      myVectCpy(memos.at(i), sMemo.obMemo(i));
    mPos = sMemo.obMpos();
  }

  void memor::adMemo(const int size) {
    assert(size >= 0 && size < INT_MAX);
    memos.clear();
    for(int i = 0; i < size; ++i)
      memos.push_back(myVectMalloc(DIM));
    mPos = memos.size();
  }

  void memor::nwMemo(const myVect* sVect) {
    assert(memos.size() > 0);
    mPos = (mPos + 1) % memos.size();
    myVectCpy(memos.at(mPos), sVect);
  }

  void memor::mkMemo(const myVect* sVect, const int indx) {
    assert(memos.size() > 0);
    int cPos = (mPos - indx) % ((signed int) memos.size());
    myVectCpy(memos.at((cPos >= 0) ? cPos : (cPos + memos.size())),
              sVect);
  }

  const myVect* memor::obMemo(const int indx) const {
    assert(memos.size() > 0);
    int cPos = (mPos - indx) % ((signed int) memos.size());
    return memos.at((cPos >= 0) ? cPos : (cPos + memos.size()));
  }

  //class DMEMO
  pmemo::pmemo() {
    mLen = 0;
    updt = false;
  }

  template <class type>
  void pmemo::adPmem(const int size, const vector<type>* objcs) {
  //void pmemo::adPmem(const int size, const vector<partk>* objcs) {
    assert(size >= 0 && size < INT_MAX);
    memos.clear();
    for(int i = 0; i < (size * (*objcs).size() * ((*objcs).size() - 1) / 2);
        ++i)
      memos.push_back(myVectMalloc(DIM));
    mLen = size;
    mPos = mLen;
    updt = true;
  }
  template void pmemo::adPmem(const int, const vector<partk>*);

  void pmemo::nwMpos() {
    assert(mLen > 0);
    mPos = (mPos + 1) % mLen;
    updt = true;
  }
 
  template <class type>
  void pmemo::nwPmem(const myVect* sVect, const vector<type>* parts,
                     const type* atom1, const type* atom2) {
  //void pmemo::nwPmem(const myVect* sVect, const vector<partk>* parts,
  //                   const partk* atom1, const partk* atom2) {
    assert(mLen > 0);
    int pos1, pos2 = obPInd (parts, atom1, atom2);
    pos1 = (mPos + 1) % mLen;
    updt = false;
    myVectCpy(memos.at(pos1 + (pos2 * mLen)), sVect);
  }
  template void pmemo::nwPmem(const myVect*, const vector<partk>*,
                              const partk*, const partk*);
 
  template <class type>
  void pmemo::mkPmem(const myVect* sVect, const vector<type>* parts,
                     const type* atom1, const type* atom2, const int indx) {
  //void pmemo::mkPmem(const myVect* sVect, const vector<partk>* parts,
  //                   const partk* atom1, const partk* atom2, const int indx) {
    assert(mLen > 0);
    int pos1, pos2 = obPInd (parts, atom1, atom2);
    pos1 = (mPos - indx) % mLen;
    pos1 = (pos1 >= 0) ? pos1 : (pos1 + mLen); 
    myVectCpy(memos.at(pos1 + (pos2 * mLen)), sVect);
  }
  template void pmemo::mkPmem(const myVect*, const vector<partk>*,
                              const partk*, const partk*, const int);

  template <class type>
  const myVect* pmemo::obPmem(const vector<type>* parts, const type* atom1,
                              const type* atom2, const int indx) const {
  //const myVect* pmemo::obPmem(const vector<partk>* parts, const partk* atom1,
  //                            const partk* atom2, const int indx) const {
    assert((mLen > 0) && updt);
    if(!updt)
      cerr << "Memory is not updated" << endl;
    int pos1, pos2 = obPInd (parts, atom1, atom2);
    pos1 = (mPos - indx) % mLen;
    pos1 = (pos1 >= 0) ? pos1 : (pos1 + mLen); 
    return memos.at(pos1 + (pos2 * mLen));
  }
  template const myVect* pmemo::obPmem(const vector<partk>*, const partk*,
                                       const partk*, const int) const;

  //class AUXMM
  //BUG: how do I know if post and velo haven't been pre-allocated
  //somewhere else
  auxMm::auxMm(const auxMm& sMemo) {
    post = myVectMalloc((*(sMemo.obPost())).size);
    velo = myVectMalloc((*(sMemo.obVelo())).size);
    myVectCpy(post, sMemo.obPost());
    myVectCpy(velo, sMemo.obVelo());
  } 

  auxMm::auxMm(const int size) {
    assert(size >= 0 && size < INT_MAX);
    post = myVectMalloc(size);
    velo = myVectMalloc(size);
  }

  auxMm::~auxMm() {
    myVectMemdel(post);
    myVectMemdel(velo);
  }

  //BUG: if the target has not been pre-allocated, this gives segfault
  const auxMm& auxMm::operator =(const auxMm& sMemo) {
    if((*post).size != (*(sMemo.obPost())).size) {
      myVectMemdel(post);
      myVectMemdel(velo);
      post = myVectMalloc((*(sMemo.obPost())).size);
      velo = myVectMalloc((*(sMemo.obPost())).size);
    }
    myVectCpy(post, sMemo.obPost());
    myVectCpy(velo, sMemo.obVelo());
  }

  void auxMm::adMemo(const int size) {
    assert(size >= 0 && size < INT_MAX);
    myVectMemdel(post);
    myVectMemdel(velo);
    post = myVectMalloc(size);
    velo = myVectMalloc(size);
  }

  void auxMm::iniTrj(istream& iput) {
    myCoor n;
    for(int i = 0; i < (*post).size; ++i) {
      iput >> n;
      myVectSet(post, i, n);
    }
  }

  void auxMm::nxtStp() {
    myVectScl(velo, STPSZE);
    myVectAdd(post, velo);
  }

  void auxMm::wrTraj(ostream& oput) const {
    for(int i = 0; i < (*post).size; ++i)
      oput << myVectGet(post, i) << '\t';
    oput << endl;
  }

  //NON-MEMBER FUNCTIONS
}
