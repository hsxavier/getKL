#include <iostream>
#include <alm.h>
#include <healpix_map.h>
#include <healpix_map_fitsio.h>
#include <arr.h>                 // For Healpix ring weights.
#include <alm_healpix_tools.h>
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
  int Nside, Scheme, lmax=300, i, j, k, l, m, L, M;
  Healpix_Map<MAP_PRECISION> NoiseMap;
  Alm<xcomplex <ALM_PRECISION> > Nlm;
  
  CovMatrix noise;
  // Testing and debugging:
  std::complex<double> z;
  int Nentries;
  
  lmax = 1;
  noise.Alloc(lmax);
  Nentries = noise.Nentries();


  z.real(1); z.imag(0); noise.set(0,0,z);
  z.real(2); z.imag(3); noise.set(0,1,z);
  z.real(4); z.imag(5); noise.set(0,2,z);
  z.real(6); z.imag(7); noise.set(0,3,z);

  z.real(2); z.imag(-3); noise.set(1,0,z);
  z.real(8); z.imag(0); noise.set(1,1,z);
  z.real(1); z.imag(2); noise.set(1,2,z);
  z.real(3); z.imag(4); noise.set(1,3,z);

  z.real(4); z.imag(-5); noise.set(2,0,z);
  z.real(1); z.imag(-2); noise.set(2,1,z);
  z.real(9); z.imag(0); noise.set(2,2,z);
  z.real(7); z.imag(0); noise.set(2,3,z);
  
  z.real(6); z.imag(-7); noise.set(3,0,z);
  z.real(3); z.imag(-4); noise.set(3,1,z);
  z.real(7); z.imag(0); noise.set(3,2,z);
  z.real(10); z.imag(0); noise.set(3,3,z);
  
  for (i=0; i<Nentries; i++) {
    for (j=0; j<Nentries; j++) {
      cout << noise(i,j) << " ";
    }
    cout<<endl;
  }
  cout<<endl;

  Announce("Solving eigensystem...");
  noise.EigenSolve();
  Announce();

  std::vector<double> Eigenvalues;
  std::vector< std::vector< std::complex<double> > > Eigenvectors;
  
  Eigenvalues.resize(Nentries);
  Eigenvectors.resize(Nentries);
  for(i=0; i<Nentries; i++) Eigenvectors[i].resize(Nentries);

  noise.Eigenvalues2(Eigenvalues);
  noise.Eigenvectors2(Eigenvectors);
  
  cout << endl << "Eigenvalues:\n";
  for (i=0;i<Nentries; i++) cout << Eigenvalues[i] << " ";
  cout << endl;
  cout << endl << "Eigenvectors:\n";
  for (i=0; i<Nentries; i++) {
    for (j=0; j<Nentries; j++) cout << Eigenvectors[i][j] << " ";
    cout<<endl;
  }
  
  return 0;

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
  Nlm.Set(lmax,lmax);
  for(l=0; l<=lmax; l++) for (m=0; m<=l; m++) { Nlm(l,m).real(0); Nlm(l,m).imag(0); }
  map2alm(NoiseMap, Nlm, RingWeights);
  Announce();
  
  // Compute angular part of the noise covariance matrix (in spherical logarithmic waves):
  
  

  return 0;
}
