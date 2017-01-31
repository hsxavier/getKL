#include <cstddef>       // For NULL pointer.
#include "CovMatrix.hpp"
#include "Utilities.hpp"
#include <vector>


/**********************************************/
/*** Initializing and terminating functions ***/
/**********************************************/


// Default Constructor:
CovMatrix::CovMatrix() {
  allocQ     = 0;
  destroyedQ = 0;
  eigenQ     = 0;
  gslCovMatrix=NULL;
  gslEigenvectors=NULL;
  gslEigenvalues=NULL;
  gslWorkspace=NULL;
  Lmax=-1;
  Nent=0;
}


// Clear all memory:
void CovMatrix::ClearAll() {
  if(gslCovMatrix   !=NULL) { gsl_matrix_complex_free(gslCovMatrix);       gslCovMatrix=NULL; }
  if(gslEigenvectors!=NULL) { gsl_matrix_complex_free(gslEigenvectors); gslEigenvectors=NULL; }
  if(gslWorkspace   !=NULL) { gsl_eigen_hermv_free(gslWorkspace);          gslWorkspace=NULL; }
  if(gslEigenvalues !=NULL) { gsl_vector_free(gslEigenvalues);           gslEigenvalues=NULL; }
  allocQ     = 0;
  destroyedQ = 0;
  eigenQ     = 0;
}


// Set matrix from sub-indices ranges:
void CovMatrix::Alloc(int lmax) {
  // Check if object was already set, and if so, clear memory first:
  if(allocQ) {
    warning("CovMatrix::Alloc: Reseting previously set matrix.");
    ClearAll();
  }
  // Allocate memory and set covariance matrix information:
  Lmax     = lmax;
  Nent = (Lmax+1)*(Lmax+1); // m goes from -l to l.
  gslCovMatrix    = gsl_matrix_complex_alloc(Nent, Nent); 
  gslEigenvectors = gsl_matrix_complex_alloc(Nent, Nent); 
  gslWorkspace    = gsl_eigen_hermv_alloc(Nent);
  gslEigenvalues  = gsl_vector_alloc(Nent);
}



/**************************************/
/*** Internal referencing functions ***/
/**************************************/


// Translate three sub-indexes (w,l,m) into a single one:
int CovMatrix::wlm2i(int w, int l, int m) const {
  // Attention: for now we are only using indexes l,m.
  return l*(l+1)+m;
}


// Translate a single index into three sub-indexes (w,l,m):
void CovMatrix::i2wlm(int i, int *w, int *l, int *m) const {
  // Attention: for now we are only using indexes l,m.
  (*l) = (int)(sqrt(i));
  (*m) = i-(*l)*(*l+1);
}



/*********************************************************************/
/*** Functions that return information about the covariance matrix ***/
/*********************************************************************/

// Returns the number the number of columns and rows (they are the same):
int CovMatrix::Nentries() const{
  return Nent;
}


/***********************************************/
/*** Functions for accessing matrix elements ***/
/***********************************************/


// Used to return values from entries in the matrix:
std::complex<double> CovMatrix::operator()(int i, int j) const {
  std::complex<double> togo;
  gsl_complex z;
  if (destroyedQ) warning("CovMatrix::operator(): Operations destroyed cov. matrix.");
  z = gsl_matrix_complex_get(gslCovMatrix, i, j);
  togo.real(GSL_REAL(z));
  togo.imag(GSL_IMAG(z));
  return togo;
}


std::complex<double> CovMatrix::operator()(int W, int L, int M, int w, int l, int m) const {
  return operator() ( wlm2i(W,L,M), wlm2i(w,l,m) );
}

// Used to assign values to entries in the matrix:
void CovMatrix::set(int i, int j, const std::complex<double> & val) {
  gsl_complex z;
  GSL_SET_COMPLEX(&z, val.real(), val.imag()); 
  gsl_matrix_complex_set(gslCovMatrix, i, j, z);
}

void CovMatrix::set(int W, int L, int M, int w, int l, int m, const std::complex<double> & val) {
  set( wlm2i(W,L,M), wlm2i(w,l,m), val );
} 


// Copies eigenvalues to 'external' vector:
void CovMatrix::Eigenvalues2(std::vector<double> & external) const {
  int i;
  if (!eigenQ) warning("CovMatrix::Eigenvalues2: Eigenvalues not computed.");
  for (i=0; i<Nent; i++) external[i] = gsl_vector_get(gslEigenvalues, i);
}
  

// Copies eigenvectors to 'external' matrix:
void CovMatrix::Eigenvectors2(std::vector< std::vector< std::complex<double> > > & external) const {
  int i, j;
  gsl_complex z;
  if (!eigenQ) warning("CovMatrix::Eigenvectors2: Eigenvectors not computed.");
  for (i=0; i<Nent; i++)
    for (j=0; j<Nent; j++) {
      z = gsl_matrix_complex_get(gslEigenvectors, j, i);
      external[i][j].real(GSL_REAL(z));
      external[i][j].imag(GSL_IMAG(z));
    }
}



/*****************************************/
/*** Functions that process the matrix ***/
/*****************************************/


// Obtains the eigenvector and eigenvalues of the matrix, destroing part of it in the process:
void CovMatrix::EigenSolve() {
  gsl_eigen_hermv(gslCovMatrix, gslEigenvalues, gslEigenvectors, gslWorkspace);
  gsl_eigen_hermv_sort(gslEigenvalues, gslEigenvectors, GSL_EIGEN_SORT_ABS_DESC);
  destroyedQ = 1;
  eigenQ     = 1;
}

