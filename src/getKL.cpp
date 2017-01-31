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
  
  lmax = 3;
  noise.Alloc(lmax);
  Nentries = noise.Nentries();

  /*
  for (i=0; i<Nentries; i++)
    for(j=0; j<Nentries; j++) {
      z.real((double)i); z.imag((double)j);
      noise.set(i,j,z);
    }
  */

  for (L=0; L<=lmax; L++) for(M=-L; M<=L; M++) {
      for (l=0; l<=lmax; l++) for(m=-l; m<=l; m++) {
	  z.real((double)(L*(L+1)+M)); z.imag((double)(l*(l+1)+m));
	  noise.set(0,L,M,0,l,m,z);
	}
    }

  for (i=0; i<Nentries; i++) {
    for(j=0; j<Nentries; j++) {
      z = noise(i,j);
      cout << z << " ";
    }
    cout<< endl;
  }
  
  /*
  for (L=0; L<=lmax; L++) for(M=-L; M<=L; M++) {
      for (l=0; l<=lmax; l++) for(m=-l; m<=l; m++) {
	  z = noise(0,L,M,0,l,m);
	  cout << z << " ";
	}
      cout << endl;
    }
  */
  return 0;

  /*gsl_matrix_complex *gmat, *evec;
  gsl_eigen_hermv_workspace *gslw;
  gsl_complex z;
  gsl_vector *eval;
  k=4;
  gmat = gsl_matrix_complex_alloc(k, k); 
  evec = gsl_matrix_complex_alloc(k, k); 
  gslw = gsl_eigen_hermv_alloc(k);
  eval = gsl_vector_alloc(k);

  /*
  GSL_SET_COMPLEX(&z, 0, 0); 
  gsl_matrix_complex_set(gmat, 0, 0, z);
  GSL_SET_COMPLEX(&z, 0, -1); 
  gsl_matrix_complex_set(gmat, 0, 1, z);
  GSL_SET_COMPLEX(&z, 0, 1); 
  gsl_matrix_complex_set(gmat, 1, 0, z);
  GSL_SET_COMPLEX(&z, 0, 0); 
  gsl_matrix_complex_set(gmat, 1, 1, z);

  GSL_SET_COMPLEX(&z, 0, 0); 
  gsl_matrix_complex_set(gmat, 0, 0, z);
  GSL_SET_COMPLEX(&z, 1, 0); 
  gsl_matrix_complex_set(gmat, 0, 1, z);
  GSL_SET_COMPLEX(&z, 1, 0); 
  gsl_matrix_complex_set(gmat, 1, 0, z);
  GSL_SET_COMPLEX(&z, 0, 0); 
  gsl_matrix_complex_set(gmat, 1, 1, z);
  
  
  GSL_SET_COMPLEX(&z, 1, 0); 
  gsl_matrix_complex_set(gmat, 0, 0, z);
  GSL_SET_COMPLEX(&z, 2, 3); 
  gsl_matrix_complex_set(gmat, 0, 1, z);
  GSL_SET_COMPLEX(&z, 4, 5); 
  gsl_matrix_complex_set(gmat, 0, 2, z);
  GSL_SET_COMPLEX(&z, 6, 7); 
  gsl_matrix_complex_set(gmat, 0, 3, z);

  GSL_SET_COMPLEX(&z, 2, -3); 
  gsl_matrix_complex_set(gmat, 1, 0, z);
  GSL_SET_COMPLEX(&z, 8, 0);
  gsl_matrix_complex_set(gmat, 1, 1, z);
  GSL_SET_COMPLEX(&z, 1, 2); 
  gsl_matrix_complex_set(gmat, 1, 2, z);
  GSL_SET_COMPLEX(&z, 3, 4); 
  gsl_matrix_complex_set(gmat, 1, 3, z);

  GSL_SET_COMPLEX(&z, 4, -5); 
  gsl_matrix_complex_set(gmat, 2, 0, z);
  GSL_SET_COMPLEX(&z, 1, -2); 
  gsl_matrix_complex_set(gmat, 2, 1, z);
  GSL_SET_COMPLEX(&z, 9, 0); 
  gsl_matrix_complex_set(gmat, 2, 2, z);
  GSL_SET_COMPLEX(&z, 7, 0); 
  gsl_matrix_complex_set(gmat, 2, 3, z);

  GSL_SET_COMPLEX(&z, 6, -7); 
  gsl_matrix_complex_set(gmat, 3, 0, z);
  GSL_SET_COMPLEX(&z, 3, -4); 
  gsl_matrix_complex_set(gmat, 3, 1, z);
  GSL_SET_COMPLEX(&z, 7, 0); 
  gsl_matrix_complex_set(gmat, 3, 2, z);
  GSL_SET_COMPLEX(&z, 10, 0); 
  gsl_matrix_complex_set(gmat, 3, 3, z);
  

  
  for (i=0; i<k; i++) {
    for (j=0; j<k; j++) {
      z = gsl_matrix_complex_get(gmat, i, j);
      printf ("%g + %gi   ",GSL_REAL(z), GSL_IMAG(z));
    }
    cout<<endl;
  }
  cout<<endl;

  gsl_eigen_hermv(gmat, eval, evec, gslw);

  for (i=0; i<k; i++) printf ("%g  ",gsl_vector_get(eval,i));
  printf("\n");
  cout<<endl;

  for (i=0; i<k; i++) {
    for (j=0; j<k; j++) {
      z = gsl_matrix_complex_get(evec, j, i);
      printf ("%g ",GSL_REAL(z));
    }
    cout<<endl;
  }
  cout<<endl;

  for (i=0; i<k; i++) {
    for (j=0; j<k; j++) {
      z = gsl_matrix_complex_get(evec, j, i);
      printf ("%g ",GSL_IMAG(z));
    }
    cout<<endl;
  }
  cout<<endl;

  
  return 0;
  */

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
