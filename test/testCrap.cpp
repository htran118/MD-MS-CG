/********************************************************************
 * Something strange is going on with passing by value
 * The correct way to do is passing by reference. How the hell do I 
 * make this basic mistake?
 ********************************************************************/

#include <iostream>

using namespace std;

void test(float& dist) {
  cout << "2:" << dist << endl;
  dist = 1.0;
  cout << "3:" << dist << endl;
}

int main() {
  float dist;
  cout << "0:" << dist << endl;
  dist = 2.0;
  cout << "1:" << dist << endl;
  test(dist);
  cout << "4:" << dist << endl;
  return 0;
}
