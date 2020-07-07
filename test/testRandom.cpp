#include <iostream>
#include <random>
#include <cmath>

using namespace std;

int NUMROL = 10000;
int NUMBIN = 100;
int BINSZE = 0.1;

int main() {
  default_random_engine ENGINE;
  normal_distribution<float> norm(0.0, 1.0);
  int rand[NUMBIN], rSqr[NUMBIN];
  float rVal, rSqT = 0;

  for(int i = 0; i < NUMBIN; ++i)
    rand[i] = 0;

  for(int i = 0; i < NUMROL; ++i) {
    rVal = norm(ENGINE);
    if(abs(rVal) < (BINSZE * NUMBIN / 2))
      ++rand[int(floor((rVal + BINSZE * NUMBIN / 2) / BINSZE))];
    rSqT += rVal * rVal;
  }

  cout << rSqT / NUMROL << endl;
  for(int i = 0; i < NUMBIN; ++i)
    //cout << i << "\t" << rand[i] << endl;

  return 0;
}
