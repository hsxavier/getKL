#include <gsl/gsl_sf_coupling.h>
#include <cmath>
#include <stdlib.h>              // For abs for int.
#include "Utilities.hpp"
#include <complex>
#include <complex.h>


// Returns the conjugate of the radial part of Zwlm: 
std::complex<double> ZetawConj(double w, double r) {
  const double sqrt2pi = 2.506628274631;
  std::complex<double> z;
  z.real(1.5);  z.imag(w);
  return exp(z*log(r))/sqrt2pi;
}


// Returns (-1)^m:
double MinusOneToPower(int m) {
  if (m%2==0) return 1.0;
  else return -1.0;
}


// Returns the integral over all sphere of: Y_l1m1 x Y_l2m2 x Y*_l3m3: 
double ThreeYlmIntegral(int l1, int m1, int l2, int m2, int l3, int m3) {
  const double FourPi = 12.56637061435917;
  int twol1, twol2, twol3;
  twol1 = 2*l1;
  twol2 = 2*l2;
  twol3 = 2*l3;
  return sqrt((double)((twol1+1)*(twol2+1)*(twol3+1))/FourPi) * MinusOneToPower(m3) * 
    gsl_sf_coupling_3j(twol1, twol2, twol3, 2*m1, 2*m2, -2*m3) *
    gsl_sf_coupling_3j(twol1, twol2, twol3,    0,     0,    0);
}


// Specify the first noise map multipole ll that contributes to the angular covariance matrix:
int FirstEll(int L, int l, int M, int m) {
  int mLimit = abs(m-M);
  int lLimit = abs(L-l);

  // If |L-l| is the limit, return it already:
  if (lLimit>=mLimit) return lLimit;
  
  // If mm = |m-M| is the limit, check if the L+l+ll=2n condition is met:
  // If l+L is even, ll must be even:
  if ((L+l)%2==0) {
    if (mLimit%2==0) return mLimit;
    else return mLimit+1;
  }
  // If l+L is odd, ll must be odd:
  else {
    if (mLimit%2==0) return mLimit+1;
    if (mLimit%2==1) return mLimit;
  }
  error("FirstEll: could not find the first ell.");
}

