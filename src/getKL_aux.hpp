#ifndef GETKL_AUX
#define GETKL_AUX 1

#include <complex>

// Returns the conjugate of the radial part of Zwlm: 
std::complex<double> ZetawConj(double w, double r);

// Returns (-1)^m:
double MinusOneToPower(int m);

// Returns the integral over all sphere of: Y_l1m1 x Y_l2m2 x Y*_l3m3: 
double ThreeYlmIntegral(int l1, int m1, int l2, int m2, int l3, int m3);

// Specify the first noise map multipole ll that contributes to the angular covariance matrix:
int FirstEll(int L, int l, int M, int m);

#endif
