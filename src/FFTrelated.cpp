#include <fftw3.h>
#include <cmath>
#include "Utilities.hpp"
#include "interpol.h"
#include <complex>
#include <vector>

// Common constants:
const double twopi = 6.283185307179586;


// Returns the logarithmic interval between r0 and r1:
double IntervalU(double r0, double r1) {
  if (r0<0 || r1<r0) error("IntervalU: The condition 0<r0<r1 must be met.");
  return log(r1)-log(r0);
}


// Computes the constant step in u=ln(r) for a given sample of size N, ranging from u0=ln(r0) to u1=ln(r1):
double DeltaU(double r0, double r1, int N) {
  return IntervalU(r0,r1)/((double)N);
}


// Computes the radial distance from the FFTlog index n:
double n2r(int n, double u0, double du) {
  // du = DeltaU(r0,r1,N);
  // u0 = log(r0);
  return exp(u0 + n*du);
}


// Computes the wavelength w from the FFTlog index q:
double q2w(int q, double U) {
  // U = IntervalU(r0, r1);
  return twopi*((double)q)/U;
}


// Compute an array of values of the tabulated function fOrig(rOrig) for the sampling points n:
void MapRadial2fftwIn(double *fftwIn, int N, double *rOrig, double *fOrig, int NOrig, double r0, double r1) {
  int n;
  double u0, du;
  u0 = log(r0);
  du = DeltaU(r0,r1,N);
  if (n2r(  0, u0, du) < r0) warning("MapRadial2fftwIn: r(n=0) < r0, will use extrapolation.");
  if (n2r(N-1, u0, du) > r1) warning("MapRadial2fftwIn: r(n=N-1) > r1, will use extrapolation.");
  for (n=0; n<N; n++) fftwIn[n] = Interpol(rOrig, NOrig, fOrig, n2r(n, u0, du));
}


// Returns exp[ i (w-W) ln(r0) ]:
std::complex<double> WwPhase(double r0, double r1, int q) {
  std::complex<double> result;
  double arg;
  arg = -twopi*((double)q)*log(r0)/IntervalU(r0, r1);
  result.real(cos(arg));
  result.imag(sin(arg));
  return result;
}


// Returns the integral of: f(r) ZW(r) Z*w(r) r^2 over r, from r0 to r1: 
void ZZr2Tranform(double *rOrig, double *fOrig, int NOrig, double r0, double r1, int N, std::vector< std::complex <double> > & result) {
  double *fftwIn;
  fftw_plan plan;
  fftw_complex *fftwOut;
  double ZZconst;
  int q, qmax, q0;
  std::complex<double> temp;

  // Prepare for Fast Fourier Transform (allocate memory and create FFTW plan):
  fftwIn  = vector<double>(0,N-1);
  fftwOut = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  plan    = fftw_plan_dft_r2c_1d(N, fftwIn, fftwOut, FFTW_ESTIMATE);
  // Sample function on FFTW sampling points n:
  MapRadial2fftwIn(fftwIn, N, rOrig, fOrig, NOrig, r0, r1);
  // Do the Fast Fourier Transform:
  fftw_execute(plan);

  // Go 
  ZZconst = IntervalU(r0,r1)/twopi/((double)N);
  qmax = N/2;
  q0   = N/2 + (N%2-1);
  result.resize(N);
  
  // Set w=0 term, which is real: 
  result[q0+0] = ZZconst * fftwOut[0][0];
  // Set w>0 terms, which are complex (actually, if N is even, q=qmax is real):
  for (q=1; q<=qmax; q++) {
    temp.real(fftwOut[q][0]);
    temp.imag(fftwOut[q][1]);
    result[q0+q] = ZZconst * WwPhase(r0,r1,q) * temp;
  }
  // Set w<0 terms, which are the conjugates of w>0:
  // If N is even, do not use -qmax which is real:
  for (q=-q0; q<0; q++) {
    temp.real( fftwOut[-q][0]);
    temp.imag(-fftwOut[-q][1]);
    result[q0+q] = ZZconst * WwPhase(r0,r1,q) * temp;
  }
  
  // Free memory used:
  fftw_destroy_plan(plan);
  fftw_free(fftwOut);
  free_vector(fftwIn, 0, N-1);
}


