#include <iostream>
#include "../libs/partk.h"
#include "../libs/const.h"

using namespace std;
using namespace MD;

myCoor GRAVITY = 9.8;

int main() {
  partk atom;
  myCoor x[DIM], v[DIM];
  int stpNum, intgrCon;
  char ctrl = 'y';
  myVect *dummy;
  dummy = myVectMalloc(DIM);
  ostream& screen = cout;

  while(ctrl == 'y') {
    cout << "Enter total number of steps:" << endl;
    cin >> stpNum;

    cout << "Enter initial position:" << endl;
    for(int i = 0; i < DIM; ++i) {
      cin >> x[i];
      myVectSet(dummy, i, x[i]);
    }
    atom.mkPost(dummy);
 
    cout << "Enter initial velocity:" << endl;
    cout << "(remember what is actually stored in velo for each method)\n";
    for(int i = 0; i < DIM; ++i) {
      cin >> v[i];
      myVectSet(dummy, i, v[i]);
    }
    atom.mkVelo(dummy);

    for(int i = 0; i < DIM - 1; ++i)
      myVectSet(dummy, i, 0);
    myVectSet(dummy, (DIM - 1), GRAVITY);
    atom.mkAccl(dummy);

    cout << "Enter integrator method: 1 for Verlet, 2 for leapfrog,"
         << "3 for velocity Verlet." << endl;
    cin >> intgrCon;

    for(int i = 0; i < stpNum; ++i) {
      if(intgrCon == 1)
      {
        atom.mkAccl(dummy);
        atom.verlet();
        atom.wrTraj(screen);
      }
      if(intgrCon == 2)
      {
        atom.mkAccl(dummy);
        atom.lpFrgA();
        atom.lpFrgB();
        atom.wrTraj(screen);
      }
      if(intgrCon == 3)
      {
        atom.vVerlA();
        atom.mkAccl(dummy);
        atom.vVerlB();
        atom.wrTraj(screen);
      }
    }

    cout << "Again?" << endl;
    cin >> ctrl;
  }
  return 0;
}
