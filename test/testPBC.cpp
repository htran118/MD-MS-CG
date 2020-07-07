#include <iostream>
#include <cmath>

using namespace std;

int main() {
  double a[2], x[2], y[2], z[2];
  char ctrl = 'y';
  cout << "Enter edge size:" << endl;
  cin >> a[0] >> a[1];
  while(ctrl == 'y') {
    cout << "Enter first coordinates:" << endl;
    cin >> x[0] >> x[1];
    cout << "Enter second coordinates:" << endl;
    cin >> y[0] >> y[1];
    for(int i = 0; i < 2; ++i)
    {
       z[i] = x[i] - y[i];
       z[i] = z[i] - a[i] * round(z[i] / a[i]);
    }
    cout << "Minimum image distance:" << endl;
    cout << z[0] << " " << z[1] << endl;
    cout << "Again?" << endl;
    cin >> ctrl;
  }
  return 0;
}
    
