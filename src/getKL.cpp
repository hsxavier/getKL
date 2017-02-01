#include <iostream>
#include <alm.h>
#include <healpix_map.h>
#include <healpix_map_fitsio.h>
#include <arr.h>                 // For Healpix ring weights.
#include <alm_healpix_tools.h>
#include <algorithm>             // for max function.
#include "getKL_aux.hpp"
#include "fitsfunctions.hpp"     // For PrepRingWeights function.
#include "Utilities.hpp"
#include "CovMatrix.hpp"

#define MAP_PRECISION double
#define ALM_PRECISION double

//#include <gsl/gsl_matrix_complex_double.h>
//#include <gsl/gsl_eigen.h>

int main() {
  using std::cout;
  using std::endl;
  std::string str1;
  Healpix_Map<MAP_PRECISION> NoiseMap;
  Alm<xcomplex <ALM_PRECISION> > Nlm;
  CovMatrix noise;
  int Nside, Scheme, lmax=2, CovN, i, j, k, l, m, L, M, ll, mm;
  std::complex<double> z, z2;
  long long1, long2;

  // Read in noise map:
  Announce("Loading noise map...");
  str1.assign("../data/boss_survey_Ns128.fits");
  read_Healpix_map_from_fits(str1, NoiseMap);
  Nside  = NoiseMap.Nside();
  Scheme = NoiseMap.Scheme();
  if (Scheme !=RING) warning("Expecting map in the RING ordering scheme.");
  Announce();

  // Get harmonic coefficients from noise map:
  arr<double> RingWeights(2*Nside);
  PrepRingWeights(1, RingWeights, Nside);
  Announce("Get harmonic coefficients from noise map...");
  Nlm.Set(2*lmax,2*lmax);
  for(l=0; l<=lmax; l++) for (m=0; m<=l; m++) { Nlm(l,m).real(0); Nlm(l,m).imag(0); }
  map2alm(NoiseMap, Nlm, RingWeights);
  Announce();
  
  // Compute angular part of the noise covariance matrix (in spherical logarithmic waves):
  Announce("Allocating memory for covariance matrix...");
  noise.Alloc(lmax);
  CovN = noise.Nentries();
  Announce();
  
  Announce("Computing covariance matrix...");
  long2=((long)CovN)*((long)CovN);
  // LOOP over Cov. Matrix elements:
#pragma omp parallel for schedule(dynamic) private(i,j,k,L,M,l,m,z,ll,z2)
  for (long1=0; long1<long2; long1++) {
    i=(int)(long1/((long)CovN));
    j=(int)(long1%((long)CovN));
    noise.i2wlm(i,&k,&L,&M);
    noise.i2wlm(j,&k,&l,&m);
    // Compute the element given by a sum over Nlm:
    z.real(0); z.imag(0);
    for (ll=FirstEll(L,l,M,m); ll<=L+l; ll=ll+2) {
      if (M>m) z2 = MinusOneToPower(M-m) * conj(Nlm(ll,M-m));
      else     z2 = Nlm(ll,m-M);
      z += z2 * ThreeYlmIntegral(ll, m-M, L, M, l, m);
    }
    // Save element to the matrix:
    noise.set(i,j,z);
  }
  Announce();
  
  
  for (i=0; i<CovN; i++) {
    for (j=0; j<CovN; j++) {
      cout << noise(i,j).imag() << " ";
    }
    cout << endl;
  }

  Announce("Solving eigensystem...");
  noise.EigenSolve();
  Announce();

  std::vector<double> Eigenvalues;
  std::vector< std::vector< std::complex<double> > > Eigenvectors;
  
  Eigenvalues.resize(CovN);
  Eigenvectors.resize(CovN);
  for(i=0; i<CovN; i++) Eigenvectors[i].resize(CovN);

  noise.Eigenvalues2(Eigenvalues);
  noise.Eigenvectors2(Eigenvectors);
  
  cout << endl << "Eigenvalues:\n";
  for (i=0;i<CovN; i++) cout << Eigenvalues[i] << " ";
  cout << endl;
  cout << endl << "Eigenvectors:\n";
  for (i=0; i<CovN; i++) {
    for (j=0; j<CovN; j++) cout << Eigenvectors[i][j] << " ";
    cout<<endl;
  }
  

  return 0;
}
