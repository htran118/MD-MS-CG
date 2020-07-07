#include <iostream>
#include <vector>

using namespace std;

typedef double myVar;

int main() {
  int vSze, pos1, pos2;
  vector <myVar> vect;
  char ctrl = 'y';
  cout << "Choose size of vector" << endl;
  cin >> vSze;
  for(int i = 0; i < vSze; ++i)
    vect.push_back(0);
  while(ctrl == 'y') {
    cout << "Choose the first position" << endl;
    cin >> pos1;
    if(pos1 >= vSze || pos1 < 0) {
      cout << "Invalid position" << endl;
      continue;
    }
    cout << "Choose the second position" << endl;
    cin >> pos2;
    if(pos2 >= vSze || pos2 < 0 || pos1 >= pos2) {
      cout << "Invalid position" << endl;
      continue;
    }
    cout << "Continue?" << endl;
    cin >> ctrl;
  }
  return 0;
}    
