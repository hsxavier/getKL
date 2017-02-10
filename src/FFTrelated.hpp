#ifndef FFTRELATED_H
#define FFTRELATED_H 1

#include <complex>
#include <vector>

double IntervalU(double r0, double r1);
double q2w(int q, double U);
void ZZr2Tranform(double *rOrig, double *fOrig, int NOrig, double r0, double r1, int N, std::vector< std::complex <double> > & result);

#endif
