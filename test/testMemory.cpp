#include <iostream>
#include <vector>
#include <time.h>
#include <stdlib.h>
#include "../libs/const.h"
#include "../libs/partk.h"
#include "../libs/memor.h"

using namespace std;
using namespace MD;

int MEMOSZ = 5;
int PARTSZ = 4;

int main() {
  srand (time(NULL));
  memor memo;
  pmemo pMem;
  vector <partk> parts;
  myVect* dumy = myVectMalloc(DIM);
  int vPos = 0, vIni = MEMOSZ + rand() % (MEMOSZ * 10);
  for(int i = 0; i < PARTSZ; ++i)
    parts.push_back(partk());
  memo.adMemo(MEMOSZ);
  pMem.adPmem(MEMOSZ, &parts);
  cout << pMem.obNmemor() << endl;
  cout << "t=" << vIni << endl;
  for(int i = 0; i < vIni; ++i) {
    myVectSet(dumy, 0, i);
    memo.nwMemo(dumy);
    vPos = 0;
    for(int j = 0; j < PARTSZ; ++j)
      for(int k = j + 1; k < PARTSZ; ++k) {
        myVectSet(dumy, 0, (i * 100 + vPos));
        ++vPos;
        pMem.nwPmem(dumy, &parts, &(parts.at(j)), &(parts.at(k)));
      }
    pMem.nwMpos();
  }
  for(int i = 0; i < MEMOSZ; ++i) {
    cout << "memo[" << i << "]=" << myVectGet(memo.obMemo(i), 0) << endl;
    for(int j = 0; j < PARTSZ; ++j)
      for(int k = j + 1; k < PARTSZ; ++k)
        cout << "partk" << j << ", partk" << k << " -> pMem[" << i << "]="
             << myVectGet(pMem.obPmem(&parts, &(parts.at(j)),
                                      &(parts.at(k)), i),
                          0) << endl;
  }

  cout <<  "*** --- ***" << endl;
  for(int i = 0; i < pMem.obNmemor(); ++i)
    cout << "pMem[" << i << "]=" << myVectGet(pMem.obMemo(i), 0) << endl;

  //for(int i = - 5 * BINSZE; i <= 5 * BINSZE; ++i)
  //  cout << i << "  " << i/BINSZE << "  " << i%BINSZE << endl;
  return 0;
}
