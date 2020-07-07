/****************************************************************************
 * Lattice generator
 * Output file is in define/coord.txt
 ***************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include "../libs/const.h"

using namespace std;
using namespace MD;

double normal(double mean, double stddev) {
  static double n2 = 0.0;
  static int n2_cached = 0;
  if(!n2_cached) {
    double x, y, r, d, n1, n2, result;
    do {
      x = 2.0 * rand() / RAND_MAX - 1;
      y = 2.0 * rand() / RAND_MAX - 1;
      r = x * x + y * y;
    } while (r == 0.0 || r > 1.0);
    d = sqrt(-2.0 * log(r) / r);
    n1 = x * d;
    n2 = y * d;
    result = n1 * stddev + mean;
    n2_cached = 1;
    return result;
  }
  else {
    n2_cached = 0;
    return n2 * stddev + mean;
  }
}

int main() {
  int m[DIM + 1], 
      dummy, wCon, rConP, rConV,
      mCon, n, cLen,
      k[DIM];
  myCoor l[DIM], a[DIM], r, scale, temp, aSze;
  char ctrl = 'n';
  FILE *pf;

  srand(time(NULL));
  m[0] = 1;
  r = 0;

  cout << "This program generates structure and coordinate file for" << endl
       << "a homogenous system. Currently support residue type LJM," << endl
       << "CPN" << endl
       << endl;

  cout << "Choose the residue type: 1 for LJM, 2 for CPN" << endl;
  cin >> mCon;
  if(mCon == 2) {
    cout << "Choose the length of the coarse-grained polyethylene" << endl;
    cin >> cLen;
    cout << "Choose the vector between beads" << endl;
    for(int i = 0; i < DIM; ++i)
      cin >> a[i];
  }

  cout << "Choose the lattice type: 1 for rectangular" << endl;
  cin >> wCon;
  if(wCon == 1) {
    cout << "Choose the lattice vector" << endl;
    for(int i = 0; i < DIM; ++i)
      cin >> l[i];
    cout << "Choose the box size" << endl;
    for(int i = 1; i <= DIM; ++i) {
      cin >> dummy;
      m[i] = m[i - 1] * dummy;
    }
  }
  else
    return 1;

  cout << "Choose to randomize position: 1 for on" << endl;
  cin >> rConP;
  if(rConP == 1) {
    cout << "Choose randomization scale (between 0 and 1)" << endl;
    cin >> scale;
    if(scale < 0 || scale > 1)
      return 1;
  }

  cout << "Choose to randomize velocity: 1 for Boltzmann distribution"
       << endl;
  cin >> rConV;
  if(rConV == 1) {
    cout << "Choose temperature" << endl;
    cin >> temp;
  }

  pf = fopen(strcFl, "w");
  if(wCon == 1) {
    fprintf(pf, "%s 12345678901234567890123456789012345678901234567890",
                ctrlIgnr.c_str());
    fprintf(pf, "123456789\n");
    fprintf(pf, "\n%s List of residues\n", ctrlIgnr.c_str());
    fprintf(pf, "%s    NAME TYPE\n", ctrlIgnr.c_str());
    if(mCon == 1) {
      for(int i = 0; i < m[DIM]; ++i)
        fprintf(pf, "%s %-5d%s\n", ctrlRsid.c_str(), (i + 1), "LJM  ");
      fprintf(pf, "%s\n", ctrlStop.c_str());
      fprintf(pf, "\n%s List of particles\n", ctrlIgnr.c_str());
      fprintf(pf, "%s    NAME RNAM TYPE\n", ctrlIgnr.c_str()); 
      for(int i = 0; i < m[DIM]; ++i) 
        fprintf(pf, "%s %-5d%-5d%s\n", ctrlAtom.c_str(), 
                                       (i + 1), (i + 1), "LJ0  ");
      fprintf(pf, "%s\n", ctrlStop.c_str());
      fprintf(pf, "\n%s List of bonds\n", ctrlIgnr.c_str());
      fprintf(pf, "%s    LEFT RIGH\n", ctrlIgnr.c_str());
      fprintf(pf, "%s\n", ctrlStop.c_str());
      fprintf(pf, "\n%s List of angles\n", ctrlIgnr.c_str());
      fprintf(pf, "%s    LEFT MIDL RIGH\n", ctrlIgnr.c_str());
      fprintf(pf, "%s\n", ctrlStop.c_str());
    }
    else if(mCon == 2) {
      for(int i = 0; i < m[DIM]; ++i)
        fprintf(pf, "%s %-5d%s\n", ctrlRsid.c_str(), (i + 1), "CPN");
      fprintf(pf, "%s\n", ctrlStop.c_str());
      fprintf(pf, "\n%s List of particles\n", ctrlIgnr.c_str());
      fprintf(pf, "%s    NAME RNAM TYPE\n", ctrlIgnr.c_str()); 
      for(int i = 0; i < m[DIM]; ++i)
        for(int j = 0; j < cLen; ++j)
          fprintf(pf, "%s %-5d%-5d%s%d   \n", ctrlAtom.c_str(),
                                              (cLen * i + j + 1),
                                              (i + 1), "C", (j + 1));
      fprintf(pf, "%s\n", ctrlStop.c_str());
      fprintf(pf, "\n%s List of bonds\n", ctrlIgnr.c_str());
      fprintf(pf, "%s    LEFT RIGH\n", ctrlIgnr.c_str());
      for(int i = 0; i < m[DIM]; ++i)
        for(int j = 0; j < cLen - 1; ++j)
          fprintf(pf, "%s %-5d%-5d\n", ctrlBond.c_str(),
                                       (cLen * i + j + 1), (cLen * i + j + 2));
      fprintf(pf, "%s\n", ctrlStop.c_str());
      fprintf(pf, "\n%s List of angles\n", ctrlIgnr.c_str());
      fprintf(pf, "%s    LEFT MIDL RIGH\n", ctrlIgnr.c_str());
      for(int i = 0; i < m[DIM]; ++i)
        for(int j = 0; j < cLen - 2; ++j)
          fprintf(pf, "%s %-5d%-5d%-5d\n", ctrlAngl.c_str(),
                                           (cLen * i + j + 1),
                                           (cLen * i + j + 2),
                                           (cLen * i + j + 3));
      fprintf(pf, "\n%s List of dihedral angles\n", ctrlIgnr.c_str());
      fprintf(pf, "%s    UPPR LEFT RIGH LOWR\n", ctrlIgnr.c_str());
      for(int i = 0; i < m[DIM]; ++i)
        for(int j = 0; j < cLen - 3; ++j)
          fprintf(pf, "%s %-5d%-5d%-5d%-5d\n", ctrlDhed.c_str(),
                                               (cLen * i + j + 1),
                                               (cLen * i + j + 2),
                                               (cLen * i + j + 3),
                                               (cLen * i + j + 4));
      fprintf(pf, "%s\n", ctrlStop.c_str());
    }
    else if(mCon == 3) {
    }
  }
  fprintf(pf, "END\n");
  fclose(pf);

  aSze = 0;
  for(int k = 0; k < DIM; ++k)
    aSze += a[k] * a[k];
  aSze = sqrt(aSze);

  pf = fopen(coorFl, "w");
  fprintf(pf, "0\n");
  if(mCon == 1) {
    if(wCon == 1) {
      for(int i = 0; i < m[DIM]; ++i) {
        dummy = i;
        for(int j = 1; j <= DIM; ++j) {
          k[DIM - j] = dummy / m[DIM - j];
          dummy -= k[DIM - j] * m[DIM - j];
        }
        for(int j = 0; j < DIM; ++j) {
          if(rConP == 1)
            r = (((double) rand() - 0.5 * RAND_MAX) / RAND_MAX * scale);
          fprintf(pf, crdFmt, (k[j] + r) * l[j]);
        }
        for(int j = 0; j < DIM; ++j)
          fprintf(pf, crdFmt, 0.0);
        fprintf(pf, "\n");
      }
    }
  }
  else if(mCon == 2) {
    if(wCon == 1) {
      for(int i = 0; i < m[DIM]; ++i) {
        dummy = i;
        for(int j = 1; j <= DIM; ++j) {
          k[DIM - j] = dummy / m[DIM - j];
          dummy -= k[DIM - j] * m[DIM - j];
        }
        for(int pos = 0; pos < cLen; ++pos) {
          r = ((double) (rand() % 1000)) / 1000;
          r = - 43.35729827*pow(r,6) + 171.8242395*pow(r,5)
              - 253.780279*pow(r,4) + 178.1903533*pow(r,3)
              - 62.3982092*pow(r,2) + 11.48466302*pow(r,1) + 0.3440359924;
          //r = 10.10448567*pow(r,6) - 12.15094332*pow(r,5)
          //  - 10.54085798*pow(r,4) + 22.35514299*pow(r,3)
          //  - 11.89902898*pow(r,2) + 3.201057093*r + 0.05345500502;
          for(int j = 0; j < DIM; ++j) {
            //if(rConP == 1)
            //  r = ((double) rand() - 0.5 * RAND_MAX) / RAND_MAX * scale;
            //else
            //  r = 0;
            //fprintf(pf, crdFmt, k[j] * l[j] + (pos + r) * a[j]);
            if(pos == 0)
              fprintf(pf, crdFmt, k[j] * l[j]);
            else {
              r /= aSze;
              fprintf(pf, crdFmt, k[j] * l[j] + pos * r * a[j]);
            }
          }
          for(int j = 0; j < DIM; ++j) {
            if(rConV == 1)
              r = sqrt(temp * BOLTZM) * normal(0.0, 1.0); 
            else
              r = 0;
            fprintf(pf, crdFmt, r);
          }
          fprintf(pf, "\n");
        }
      }
    }
  }
  else if(mCon == 3) {
  }
  fclose(pf);

  return 0;
}   
