#ifndef COVMATRIX_H
#define COVMATRIX_H 1

#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_eigen.h>
#include <complex>
#include <vector>

// Indicators for what data to output in Print function:
enum DataOut {CovOut, EvalOut, EvecOut};

// Class declaration:
class CovMatrix {
private:
  bool allocQ, destroyedQ, eigenQ;
  gsl_matrix_complex *gslCovMatrix, *gslEigenvectors;
  gsl_eigen_hermv_workspace *gslWorkspace;
  gsl_vector *gslEigenvalues;
  int Lmax;
  int Nent;
public:
  CovMatrix();
  // Set sizes and allocate memory:
  void Alloc(int lmax);
  // Returns the size N of the matrix NxN:
  int Nentries() const;
  // Used to return values from entries in the matrix:
  std::complex<double> operator()(int i, int j) const;
  std::complex<double> operator()(int W, int L, int M, int w, int l, int m) const;
  // Used to assign values to entries in the matrix:
  void set(int i, int j, const std::complex<double> & val);
  void set(int W, int L, int M, int w, int l, int m, const std::complex<double> & val);
  // Obtains the eigenvalues and eigenvectors:
  void EigenSolve();
  // Returns eigenvalues to a std::vector:
  void Eigenvalues2(std::vector<double> & external) const;
  // Returns the ith eigenvector:
  void Eigenvectors2(std::vector< std::vector< std::complex<double> > > & external) const;
  // Deallocate memory:
  void ClearAll();
  // Translate sub-indices (w,l,m) into a single covariance matrix index i:
  int  wlm2i(int  w, int  l, int  m) const;
  void i2wlm(int i, int *w, int *l, int *m) const;
  // Print data to file:
  void Print(std::string fileprefix, DataOut data);
};
  
#endif
  
