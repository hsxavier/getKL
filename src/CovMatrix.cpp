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
  Lmax = -1;
  Qmax = -1;
  Nang =  0;
  Nent =  0;
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


// Set matrix from sub-indices ranges, using only spherical harmonic part:
void CovMatrix::Alloc(int lmax) {
  // Check if object was already set, and if so, clear memory first:
  if(allocQ) {
    warning("CovMatrix::Alloc: Reseting previously set matrix.");
    ClearAll();
  }
  // Allocate memory and set covariance matrix information:
  Lmax = lmax;
  Nent = (Lmax+1)*(Lmax+1); // m goes from -l to l.
  gslCovMatrix    = gsl_matrix_complex_alloc(Nent, Nent); 
  gslEigenvectors = gsl_matrix_complex_alloc(Nent, Nent); 
  gslWorkspace    = gsl_eigen_hermv_alloc(Nent);
  gslEigenvalues  = gsl_vector_alloc(Nent);
}


// Set matrix from sub-indices ranges:
void CovMatrix::Alloc(int qmax, int lmax) {
  // Check if object was already set, and if so, clear memory first:
  if(allocQ) {
    warning("CovMatrix::Alloc: Reseting previously set matrix.");
    ClearAll();
  }
  // Allocate memory and set covariance matrix information:
  Lmax = lmax;
  Nang = (Lmax+1)*(Lmax+1);
  Qmax = qmax;
  Nent = (2*Qmax+1)*Nang; // Nq=(2Qmax+1) blocks for which m goes from -l to l.
  gslCovMatrix    = gsl_matrix_complex_alloc(Nent, Nent); 
  gslEigenvectors = gsl_matrix_complex_alloc(Nent, Nent); 
  gslWorkspace    = gsl_eigen_hermv_alloc(Nent);
  gslEigenvalues  = gsl_vector_alloc(Nent);
}




/**************************************/
/*** Internal referencing functions ***/
/**************************************/


// Translate three sub-indexes (q,l,m) into a single one:
int CovMatrix::qlm2i(int q, int l, int m) const {
  if (q < -Qmax) warning("qlm2i: q < -Qmax");
  if (q >  Qmax) warning("qlm2i: q > Qmax");
  if (l <  0   ) warning("qlm2i: l < 0");
  if (l >  Lmax) warning("qlm2i: l > Lmax");
  if (m < -l   ) warning("qlm2i: m < -l");
  if (m >  l   ) warning("qlm2i: m > l");
  return (q+Qmax)*Nang + l*(l+1)+m;
}


// Translate a single index into three sub-indexes (q,l,m):
void CovMatrix::i2qlm(int i, int *q, int *l, int *m) const {
  int j;
  if (i<0)     warning("i2qlm: i < 0");
  if (i>=Nent) warning("i2qlm: i > # of entries");
  *q = i/Nang - Qmax;
  j  = i%Nang;
  (*l) = (int)(sqrt(j));
  (*m) = j-(*l)*(*l+1);
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


std::complex<double> CovMatrix::operator()(int Q, int L, int M, int q, int l, int m) const {
  return operator() ( qlm2i(Q,L,M), qlm2i(q,l,m) );
}

// Used to assign values to entries in the matrix:
void CovMatrix::set(int i, int j, const std::complex<double> & val) {
  gsl_complex z;
  GSL_SET_COMPLEX(&z, val.real(), val.imag()); 
  gsl_matrix_complex_set(gslCovMatrix, i, j, z);
}

void CovMatrix::set(int Q, int L, int M, int q, int l, int m, const std::complex<double> & val) {
  set( qlm2i(Q,L,M), qlm2i(q,l,m), val );
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



/********************************/
/*** Functions to output data ***/
/********************************/


// Print data to file:
void CovMatrix::Print(std::string fileprefix, DataOut data) {
  using std::endl;
  std::ofstream outfile1, outfile2; 
  gsl_matrix_complex *gslMat;
  int i, j;

  // Real vector case:
  if (data==EvalOut) {
    // Check if can open file for output:
    outfile1.open(fileprefix.c_str());
    if (!outfile1.is_open()) warning("CovMatrix::Print: cannot open file "+fileprefix);
    else {
      // Write eigenvalues to file:
      if (!eigenQ) warning("CovMatrix::Print: Eigenvalues not computed.");
      for (i=0; i<Nent; i++) outfile1 << gsl_vector_get(gslEigenvalues, i) << endl; 
      outfile1.close();
    }
  }

  // Complex matrix case:
  if (data==EvecOut || data==CovOut) {
    // Check if can open file for output:
    outfile1.open((fileprefix+"-re").c_str());
    if (!outfile1.is_open()) warning("CovMatrix::Print: cannot open file "+fileprefix+"-re");
    outfile2.open((fileprefix+"-im").c_str());
    if (!outfile2.is_open()) warning("CovMatrix::Print: cannot open file "+fileprefix+"-im");
    else {
      // Set output for Eigenvectors or Covariance matrix:
      if (data==EvecOut) {
	gslMat=gslEigenvectors;
	if (!eigenQ) warning("CovMatrix::Print: Eigenvectors not computed.");
      }
      else if (data==CovOut) {
	gslMat=gslCovMatrix;
	if (destroyedQ) warning("CovMatrix::Print: Covariance matrix was destroyed.");
      }
      // Write matrix to files (real and imag parts in different files):
      for (i=0; i<Nent; i++) {
	for (j=0; j<Nent; j++) {
	  outfile1 << GSL_REAL(gsl_matrix_complex_get(gslMat, i, j)) << " ";
	  outfile2 << GSL_IMAG(gsl_matrix_complex_get(gslMat, i, j)) << " ";
	}
	outfile1 << endl;
	outfile2 << endl;
      }
      outfile1.close();
      outfile2.close();
    }
  }
}
