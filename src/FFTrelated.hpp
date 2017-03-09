#ifndef FFTRELATED_H
#define FFTRELATED_H 1

#include <complex>
#include <vector>
#include "ParameterList.hpp"

double IntervalU(double r0, double r1);
double q2w(int q, double U);
int qq2i(int q, int N);
void ZZr2Tranform(double *rOrig, double *fOrig, int NOrig, double r0, double r1, 
		  double rStart, double rEnd, int N, std::vector< std::complex <double> > & result);
void PrintZZr2(std::string filename, double rStart, double rEnd, int Nfft, const std::vector<std::complex<double> > & trafo);
#endif
