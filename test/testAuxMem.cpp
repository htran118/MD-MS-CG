#include <iostream>
#include <vector>
#include <gsl/gsl_vector_float.h>
//#include "../libs/memor.h"
//#include "../libs/const.h"

//using namespace MD;
using namespace std;

class auxMm
{
  public:
    auxMm(const int size) {
      for(int i = 0; i < size; ++i) {
        post = gsl_vector_float_calloc(10);
        velo = gsl_vector_float_calloc(10);
      }
    }
  private:
    gsl_vector_float *post, *velo;
};

int main() {
  vector<auxMm> axmem;
  for(int i = 0; i < 10; ++i) {
    axmem.push_back(auxMm(999999));
    cout << "FML" << endl;
  }
  return 0;
}
