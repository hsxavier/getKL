#include <fftw3.h>
#include <cmath>
#include "Utilities.hpp"
#include "interpol.h"
#include <complex>
#include <vector>
#include "ParameterList.hpp"
#include "Integral.hpp"

/*** NOTE: Here, qq = (Q-q) is the difference between two radial mode indices Q and q. ***/

// Common constants:
const double twopi = 6.283185307179586;



/**********************************/
/*** General acessory functions ***/
/**********************************/

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
// The offset goes 0<=offset<1 and is the position inside the bin of width du.
double n2r(int n, double uStart, double du, double offset) {
  // du = DeltaU(rStart,rEnd,N);
  // uStart = log(rStart);
  return exp(uStart + (n+offset)*du);
}


// Computes the wavelength w from the FFTlog index q:
double q2w(int q, double U) {
  // U = IntervalU(rStart, rEnd);
  return twopi*((double)q)/U;
}


// Computes the position in a vector of FFTlog index qq=Q-q that goes from -qqmax to +qqmax:
int qq2i(int qq, int Nfft) {
  int qqmax;
  
  qqmax = Nfft/2;
  if (qq<-qqmax) warning("qq2i: Q-q is less than expected.");
  if (qq> qqmax) warning("qq2i: Q-q is more than expected.");
  return qq + qqmax;
}


// Print transform to file:
void PrintZZr2(std::string filename, double rStart, double rEnd, int Nfft, const std::vector<std::complex<double> > & trafo) {
  using std::endl;
  int qq, qqmax;
  double du;
  std::ofstream outfile;
  
  // Set variables:
  qqmax = Nfft/2;
  du    = IntervalU(rStart,rEnd);
  if (trafo.size() != 2*qqmax+1) warning("PrintZZr2: transform size is not in harmony with FFT vector size.");

  // Open file:
  outfile.open(filename.c_str());
  if (!outfile.is_open()) warning("PrintZZr2: cannot open file "+filename);
  else {
    // Write trafo table:
    outfile << "# Q-q, W-w, [Re], [Im]" << endl;  
    for (qq=-qqmax; qq<=qqmax; qq++) {
      outfile<< qq <<" "<< q2w(qq,du) <<" "<< trafo[qq2i(qq,Nfft)].real() <<" "<< trafo[qq2i(qq,Nfft)].imag() << endl;
    }
    outfile.close();
  }  
}


/**************************************/
/*** Functions for the FFT approach ***/
/**************************************/

// Compute an array of values of the tabulated function fOrig(rOrig) for the sampling points n:
void MapRadial2fftwIn(double *fftwIn, int N, double *rOrig, double *fOrig, int NOrig, double r0, double r1, 
		      double rStart, double rEnd) {
  const double rEPS = 1e-8;
  int n;
  double uStart, du, r;
  uStart = log(rStart);
  du = DeltaU(rStart,rEnd,N);

  if (r0<rOrig[0]) warning("MapRadial2fftwIn: requested non-zero interval starts before tabulated values, will use extrapolation.");
  if (r1>rOrig[NOrig-1]) warning("MapRadial2fftwIn: requested non-zero interval starts before tabulated values, will use extrapolation.");
  
  for (n=0; n<N; n++) {
    r = n2r(n,uStart,du,0);
    if (r<r0*(1.0-rEPS) || r>r1*(1.0+rEPS)) fftwIn[n] = 0.0; // Apply radial window to f(r).
    else fftwIn[n] = Interpol(rOrig, NOrig, fOrig, r);       // Interpolate f(r) to sampled value.
  }
}


// Returns exp[ i (w-W) ln(rStart) ]:
std::complex<double> WwPhase(double rStart, double rEnd, int qq) {
  std::complex<double> result;
  double arg;
  arg = -twopi*((double)qq)*log(rStart)/IntervalU(rStart, rEnd);
  result.real(cos(arg));
  result.imag(sin(arg));
  return result;
}


// Returns the integral of: f(r) Z_W(r) Z*_w(r) r^2 over r, from rStart to rEnd, but zeroing f(r) outside [r0,r1]:
// It returns the results for the each W-w.    !! USING FFT !!
void ZZr2TranformFFT(double *rOrig, double *fOrig, int NOrig, double r0, double r1, 
		     double rStart, double rEnd, int Nfft, std::vector< std::complex <double> > & result) {
  double *fftwIn;
  fftw_plan plan;
  fftw_complex *fftwOut;
  double ZZconst;
  int qq, qqmax;             // qq = Q-q.
  std::complex<double> temp;

  // Set result range:
  qqmax = Nfft/2;                 // This is the max. qq used to compute the covariance, which is twice the max. q for the radial modes.
  result.resize(Nfft + 1-Nfft%2); // We will make the result cover a simmetric range, i.e., go from -qqmax to qqmax.
                                  // Therefore, it will always be odd-sized.

  // Prepare for Fast Fourier Transform (allocate memory and create FFTW plan):
  fftwIn  = vector<double>(0,Nfft-1);
  fftwOut = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nfft);
  plan    = fftw_plan_dft_r2c_1d(Nfft, fftwIn, fftwOut, FFTW_ESTIMATE);

  // Sample function on FFTW sampling points n:
  MapRadial2fftwIn(fftwIn, Nfft, rOrig, fOrig, NOrig, r0, r1, rStart, rEnd);
  // Do the Fast Fourier Transform:
  fftw_execute(plan);

  // Combine Integral with other factors:
  ZZconst = DeltaU(rStart,rEnd,Nfft)/twopi;
  // Set W-w=0 term, which is real: 
  result[qq2i(0,Nfft)] = ZZconst * fftwOut[0][0];
  // Set W-w>0 terms, which are complex (actually, if Nfft is even, fftwOut with qq=qqmax is real):
  for (qq=1; qq<=qqmax; qq++) {
    temp.real(fftwOut[qq][0]);
    temp.imag(fftwOut[qq][1]);
    result[qq2i(qq,Nfft)] = ZZconst * WwPhase(rStart,rEnd,qq) * temp;
  }
  // Set W-w<0 terms, which are the conjugates of W-w>0:
  for (qq=-qqmax; qq<0; qq++) {
    temp.real( fftwOut[-qq][0]);
    temp.imag(-fftwOut[-qq][1]);
    result[qq2i(qq,Nfft)] = ZZconst * WwPhase(rStart,rEnd,qq) * temp;
  }
  
  // Free memory used:
  fftw_destroy_plan(plan);
  fftw_free(fftwOut);
  free_vector(fftwIn, 0, Nfft-1);
}



/*****************************************/
/*** Functions for the Romberg aproach ***/
/*****************************************/

// Integrand of the radial part of the noise covariance matrix (ZZr2 transform):
// (it returns either real or imag. part if iPhase=0 or -pi/2, respectively).
double ZZr2Integrand(double iPhase, double *rArray, double *noiseArray, int Na, double wdiff, double r) {
  return Interpol(rArray, Na, noiseArray, r) * cos(-wdiff*log(r)+iPhase) / (2*M_PI) / r;
}
double ZZr2IntegrandRe(double *rArray, double *noiseArray, int Na, double wdiff, double r) {
  return ZZr2Integrand(0, rArray, noiseArray, Na, wdiff, r);
}
double ZZr2IntegrandIm(double *rArray, double *noiseArray, int Na, double wdiff, double r) {
  return ZZr2Integrand(-M_PI_2, rArray, noiseArray, Na, wdiff, r);
}


// Solve ZZr2 noise integral for a single W-w 'wdiff'.
// rArray[0...Na-1] and noiseArray[0...Na-1] are tabulated density contrast noise [1/n(r)];
// The integral is performed from r0 to r1, given by the radial selection function limits.
void ZZr2Romberg(double *rArray, double *noiseArray, int Na, double wdiff, double r0, double r1, 
		 double *ResRe, double *ResIm) {
  *ResRe = qrombTable(ZZr2IntegrandRe, r0, r1, rArray, noiseArray, Na, wdiff); 
  *ResIm = qrombTable(ZZr2IntegrandIm, r0, r1, rArray, noiseArray, Na, wdiff);
}


// Returns the integral of: f(r) Z_W(r) Z*_w(r) r^2 over r, from rStart to rEnd, but zeroing f(r) outside [r0,r1]:
// It returns the results for the each W-w.    !! USING Romberg !!
void ZZr2TranformRomb(double *rOrig, double *fOrig, int NOrig, double r0, double r1, 
		      double rStart, double rEnd, int Nfft, std::vector< std::complex <double> > & result) {
  int qq, qqmax;
  double du, ResRe, ResIm;

  qqmax = Nfft/2;                 // Check ZZr2Tranform for more information on this operation.
  du    = IntervalU(rStart,rEnd);

  result.resize(Nfft + 1-Nfft%2); // Check ZZr2Tranform for more information on this operation.

  //#pragma omp parallel for private(ResRe, ResIm)
  for (qq=-qqmax; qq<=qqmax; qq++) {
    ZZr2Romberg(rOrig, fOrig, NOrig, q2w(qq,du), r0, r1, &ResRe, &ResIm);
    result[qq2i(qq,Nfft)].real(ResRe);
    result[qq2i(qq,Nfft)].imag(ResIm);
  }

}
