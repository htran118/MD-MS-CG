/****************************************************************************
 * Constants
 ***************************************************************************/

#ifndef CONST_H
#define CONST_H

#include <string>
#include <random>
#include <chrono>
#include <iostream>
#include <gsl/gsl_vector_float.h>
#include <gsl/gsl_blas.h>

using namespace std;

namespace MD
{
  //GLOBAL ALIASES
  typedef float myCoor;
  typedef gsl_vector_float myVect;
  typedef gsl_matrix_float myMatr;

  static const auto& myVectMalloc = gsl_vector_float_calloc;
  static const auto& myVectMemdel = gsl_vector_float_free;
  static const auto& myVectSet = gsl_vector_float_set;
  static const auto& myVectGet = gsl_vector_float_get;
  static const auto& myVectMin = gsl_vector_float_min;
  static const auto& myVectMax = gsl_vector_float_max;
  static const auto& myVectCpy = gsl_vector_float_memcpy;
  static const auto& myVectZro = gsl_vector_float_set_zero;

  static const auto& myMatrMalloc = gsl_matrix_float_calloc;
  static const auto& myMatrMemdel = gsl_matrix_float_free;
  static const auto& myMatrSet = gsl_matrix_float_set;
  static const auto& myMatrSetrow = gsl_matrix_float_set_row;
  static const auto& myMatrGet = gsl_matrix_float_get;
  static const auto& myMatrZro = gsl_matrix_float_set_zero;

  static const auto& myVectAdd = gsl_vector_float_add;
  static const auto& myVectSub = gsl_vector_float_sub;
  static const auto& myVectScl = gsl_vector_float_scale;
  static const auto& myVectMul = gsl_vector_float_mul;
  static const auto& myVectDiv = gsl_vector_float_div;
  static const auto& myMatrAdd = gsl_matrix_float_add;
  static const auto& myMatrScl = gsl_matrix_float_scale;
  static const auto& myVectDot = gsl_blas_sdot;
  static const auto& myVectNrm = gsl_blas_snrm2;
  static const auto& myMatrVecmul = gsl_blas_sgemv; 
  //GLOBAL VARIABLES
  static ostream& COUTPT = cout;
  static istream& CINPUT = cin;
  extern const myCoor ONEPIE, TWOPIE, ROOTWO, IVRTWO, COULMB, BOLTZM, MINZRO;
  extern const myCoor STPSZE, STPHLF, STPINV, STPIHF, STEPSQ, STPIRT;
  extern const int DIM;
  extern const int MAXNPART, MAXNBOND, MAXNANGL, MAXNDHED;
  extern const int BONDFORC;
  extern const bool ANGLPART, LJONPART, DPDPARTK, LONGMEMO;
  extern const char *paraFl, *topoFl, *strcFl, *coorFl, *trajFl, *memoFl,
                    *convFl, *defNme, *datNme, *cgStFl;
  extern const char *crdFmt, *typFmt, *nmeFmt, *sNmFmt;
  extern const string ctrlAtom, ctrlMass, ctrlBond, ctrlAngl, ctrlDhed, 
                      ctrlRsid, ctrlTime, ctrlStop, ctrlIgnr, ctrlEnd,
                      emptyNme;
  //GLOBAL FUNCTIONS
  static default_random_engine ENGINE;
  static normal_distribution<myCoor> NRMDIS{0.0, 1.0};
  static uniform_real_distribution<myCoor> UNIDIS(0.0,1.0);
  void myVectInv(myVect* vec1, const myCoor leng);
  void myVectRev(myVect* vec1);
  void myVectCrs(myVect* vec1, const myVect* vec2);
  void myVectVecCrs(const myVect* vec1, const myVect* vec2, myVect* vec3);
  void myVectUnt(myVect* uVec, myCoor& dist);
  void myVectRnd(myVect* uVec);
  void myVectOrd(myVect* sVec, const myVect* uVec);
  void myVectInp(myVect* vect, istream& iput);
  void rndVec(myVect* uVec);
  void myVectOut(const myVect* vect);
  void dbuggr();
  //GLOBAL DUMMY VECTOR
  static myVect *FOR1 = myVectMalloc(DIM), *FOR2 = myVectMalloc(DIM),
                *FOR3 = myVectMalloc(DIM), *FTOT = myVectMalloc(DIM),
                *FORC = myVectMalloc(DIM), 
                *VELO = myVectMalloc(DIM), *POST = myVectMalloc(DIM),
                *PBC1 = myVectMalloc(DIM), *PBC2 = myVectMalloc(DIM);
  //NOTE: be careful when using these vectors, as they may be 
  //      used elsewhere
  static myVect *NRM1 = myVectMalloc(DIM), *NRM2 = myVectMalloc(DIM),
                *LEN1 = myVectMalloc(DIM), *LEN2 = myVectMalloc(DIM),
                *LEN3 = myVectMalloc(DIM), *UNIT = myVectMalloc(DIM);
}
#endif
