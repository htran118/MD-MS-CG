#include <iostream>

using namespace std;

int main() {
  unsigned int s = 5;
  auto d = (signed int) s;
  cout << "d=" << d << endl;
  for(int i = -15; i < 15; ++i)
    cout << i << "%s=" << i%s << "," << i << "%d=" << i%d <<  endl;
  return 0;
}
