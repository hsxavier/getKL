#include <fftw3.h>
#include <cmath>
#include "Utilities.hpp"
#include "interpol.h"
#include <complex>
#include <vector>
#include "ParameterList.hpp"

// Common constants:
const double twopi = 6.283185307179586;


// Returns the logarithmic interval between rStart and rEnd:
double IntervalU(double rStart, double rEnd) {
  if (rStart<0 || rEnd<rStart) error("IntervalU: The condition 0<rStart<rEnd must be met.");
  return log(rEnd)-log(rStart);
}


// Computes the constant step in u=ln(r) for a given sample of size N, ranging from uStart=ln(rStart) to uEnd=ln(rEnd):
double DeltaU(double rStart, double rEnd, int N) {
  return IntervalU(rStart,rEnd)/((double)N);
}


// Computes the radial distance from the FFTlog index n:
double n2r(int n, double uStart, double du) {
  // du = DeltaU(rStart,rEnd,N);
  // uStart = log(rStart);
  return exp(uStart + n*du);
}


// Computes the wavelength w from the FFTlog index q:
double q2w(int q, double U) {
  // U = IntervalU(rStart, rEnd);
  return twopi*((double)q)/U;
}


// Compute an array of values of the tabulated function fOrig(rOrig) for the sampling points n:
void MapRadial2fftwIn(double *fftwIn, int N, double *rOrig, double *fOrig, int NOrig, double r0, double r1, 
		      double rStart, double rEnd) {
  int n;
  double uStart, du, r;
  uStart = log(rStart);
  du = DeltaU(rStart,rEnd,N);

  if (r0<rOrig[0]) warning("MapRadial2fftwIn: requested non-zero interval starts before tabulated values, will use extrapolation.");
  if (r1>rOrig[NOrig-1]) warning("MapRadial2fftwIn: requested non-zero interval starts before tabulated values, will use extrapolation.");
  
  for (n=0; n<N; n++) {
    r = n2r(n,uStart,du);
    if (r<r0 || r>r1) fftwIn[n] = 0.0;                  // Apply radial window to f(r).
    else fftwIn[n] = Interpol(rOrig, NOrig, fOrig, r);  // Interpolate f(r) to sampled value.
  }
}


// Returns exp[ i (w-W) ln(rStart) ]:
std::complex<double> WwPhase(double rStart, double rEnd, int q) {
  std::complex<double> result;
  double arg;
  arg = -twopi*((double)q)*log(rStart)/IntervalU(rStart, rEnd);
  result.real(cos(arg));
  result.imag(sin(arg));
  return result;
}


// Returns the integral of: f(r) ZW(r) Z*w(r) r^2 over r, from rStart to rEnd, but zeroing f(r) outside [r0,r1]: 
void ZZr2Tranform(double *rOrig, double *fOrig, int NOrig, double r0, double r1, 
		  double rStart, double rEnd, int N, std::vector< std::complex <double> > & result) {
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
  MapRadial2fftwIn(fftwIn, N, rOrig, fOrig, NOrig, r0, r1, rStart, rEnd);

  // debug:
  int i;
  printf("fftwIn:\n");
  for(i=0; i<N; i++) printf("%g ",fftwIn[i]);
  printf("\n");

  // Do the Fast Fourier Transform:
  fftw_execute(plan);

  // debug: 
  printf("fftwOut[Re]:\n");
  for(i=0; i<N/2+1; i++) printf("%g ",fftwOut[i][0]);
  printf("\n");
  printf("fftwOut[Im]:\n");
  for(i=0; i<N/2+1; i++) printf("%g ",fftwOut[i][1]);
  printf("\n");



  // Go 
  ZZconst = DeltaU(rStart,rEnd,N)/twopi;
  qmax = N/2;
  q0   = N/2 + (N%2-1);
  result.resize(N);
  
  // Set w=0 term, which is real: 
  result[q0+0] = ZZconst * fftwOut[0][0];
  // Set w>0 terms, which are complex (actually, if N is even, q=qmax is real):
  for (q=1; q<=qmax; q++) {
    temp.real(fftwOut[q][0]);
    temp.imag(fftwOut[q][1]);
    result[q0+q] = ZZconst * WwPhase(rStart,rEnd,q) * temp;
  }
  // Set w<0 terms, which are the conjugates of w>0:
  // If N is even, do not use -qmax which is real:
  for (q=-q0; q<0; q++) {
    temp.real( fftwOut[-q][0]);
    temp.imag(-fftwOut[-q][1]);
    result[q0+q] = ZZconst * WwPhase(rStart,rEnd,q) * temp;
  }
  
  // Free memory used:
  fftw_destroy_plan(plan);
  fftw_free(fftwOut);
  free_vector(fftwIn, 0, N-1);
}


