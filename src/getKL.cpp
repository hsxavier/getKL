/* Code for computing the Karhunen-Loeve (KL) modes for a galaxy survey 
   with separable selection function. The idea is that we first compute 
   spherical logarithmic waves modes (angular part are expanded in 
   spherical harmonics and the radial part are logarithmic waves), and 
   then we find the eigenvalues and eigenvectors of the covariance matrix 
   to build a reduced set of linear combinations of these modes that best 
   describe the data, thus obtaining the KL modes. 
*/

#include <iostream>
#include <alm.h>
#include <healpix_map.h>
#include <healpix_map_fitsio.h>
#include <arr.h>                 // For Healpix ring weights.
#include <alm_healpix_tools.h>   // For map2alm.
#include "getKL_aux.hpp"
#include "fitsfunctions.hpp"     // For PrepRingWeights function.
#include "Utilities.hpp"
#include "CovMatrix.hpp"
#include "ParameterList.hpp"
#include "definitions.h"
#include "output.hpp"



/*********************************************************************************/
/***                                 Main code                                 ***/
/*** Compute covariance matrix for harmonic coefficients of a Map              ***/
/*** In other words:                                                           ***/
/***          J_lmLM = \sum_l'm' Al'm' \int Yl'm'(r) Ylm(r) Y*LM(r) d2r        ***/  
/*********************************************************************************/

int main(int argc, char *argv[]) {
  using std::cout;
  using std::endl;
  using namespace ParDef; ParameterList config;         // Easy configuration file use.
  // Special variables:
  time_t StartAll;                                      // For timing the code run.
  Healpix_Map<MAP_PRECISION> Map;
  Alm<xcomplex <ALM_PRECISION> > Alm;
  CovMatrix angular;
  int Nside, Scheme, lmax, AngN;
  // Generic variables and iterators:
  std::string str1;
  int i, j, q, Q, l, m, L, M, ll, mm;
  std::complex<double> z, z2;
  long long1, long2;

  StartAll = time(NULL);


  /***************************************/
  /*** Part 0: Load configuration file ***/
  /***************************************/

  cout << endl;
  cout << "Code commit:  " << COMMIT << endl;
  // Loading config file:
  if (argc<=1) { cout << "You must supply a config file." << endl<<endl; return 0;}
  config.load(argv[1]);
  cout << endl;
  cout << "-- Configuration setup:\n";
  cout << "   File: "<<argv[1]<<endl;
  config.lineload(argc, argv);
  config.show();
  cout << endl;
  lmax = config.readi("LMAX");
  

  /********************************************************/
  /*** Part 1: Load and process selection function data ***/
  /********************************************************/

  /*** Part 1.2: Prepare angular data ***/
  cout << "** Prepare angular data:"<<endl;

  // Read in completeness map and compute density contrast noise map:
  Announce("   Loading angular map...");
  read_Healpix_map_from_fits(config.reads("ANGULAR_MAP"), Map);
  Nside  = Map.Nside();
  Scheme = Map.Scheme();
  if (Scheme !=RING) warning("Expecting map in the RING ordering scheme.");
  Announce();
  cout << "   ANGULAR_MAP has Nside = "<<Nside<<".\n";

  // Get harmonic coefficients from noise map:
  arr<double> RingWeights(2*Nside);
  PrepRingWeights(1, RingWeights, Nside);
  Announce("   Get harmonic coeffs up to 2LMAX from noise map...");
  Alm.Set(2*lmax,2*lmax); // Noise cov. matrix depends on integral of 3 Ylm up to Alm with l=2*LMAX.
  for(l=0; l<=2*lmax; l++) for (m=0; m<=l; m++) { Alm(l,m).real(0); Alm(l,m).imag(0); }
  map2alm(Map, Alm, RingWeights);
  Announce();
  WriteAlm(Alm, config, "ALM_OUT", 1);


  /***********************************************/
  /*** Part 2: Compute noise covariance matrix ***/
  /***********************************************/


  /*** Part 2.1: Compute angular part of the noise covariance matrix ***/
  cout << "** Prepare angular covariance matrix:"<<endl;

  // Prepare Covariance matrix object:
  Announce("   Allocating memory for matrix...");
  angular.Alloc(lmax);
  AngN = angular.Nentries();
  Announce();

  // Compute the noise angular covariance matrix: 
  cout <<  "   Angular Cov. matrix is N="<<AngN<<endl;  
  Announce("   Computing diag. and below of ang. cov. matrix...");
  long2=(((long)AngN+1)*((long)AngN))/2;
  // LOOP over Cov. Matrix elements:
#pragma omp parallel for schedule(dynamic) private(i,j,Q,L,M,q,l,m,z,ll,z2)
  for (long1=0; long1<long2; long1++) {
    i = (int)((sqrt(8*long1+1)-1.0)/2.0);
    j = (int)(long1-(i*(i+1))/2);
    angular.i2qlm(i,&q,&l,&m);
    angular.i2qlm(j,&Q,&L,&M);
    // Compute the element given by a sum over Alm:
    z.real(0); z.imag(0);
    for (ll=FirstEll(L,l,M,m); ll<=L+l; ll=ll+2) {
      if (M>m) z2 = MinusOneToPower(M-m) * conj(Alm(ll,M-m));
      else     z2 = Alm(ll,m-M);
      z += z2 * ThreeYlmIntegral(ll, m-M, L, M, l, m);
    }
    // Save element to the matrix:
    angular.set(i,j,z);
  }
  Announce();

  // If requested, set remaining elements of cov. matrix as conjugate transpose:
  // NOTE: we are currently resetting the diagonal. Improve this in the future !!
  if (config.readi("FULLCOV")==1) {
    Announce("   Filling upper part of ang. cov. matrix...");
#pragma omp parallel for schedule(dynamic) private(i,j)
    for (long1=0; long1<long2; long1++) {
      i = (int)((sqrt(8*long1+1)-1.0)/2.0);
      j = (int)(long1-(i*(i+1))/2);
      angular.set(j,i, conj(angular(i,j)));
    }
    Announce();
  }
  
  // If requested, print out the covariance matrix:
  str1=config.reads("ANGCOV_OUT");
  if (str1!="0") {
    Announce(">> Writting angular part of Cov. Matrix to "+str1+"...");
    angular.Print(str1, CovOut);
    Announce();
  }

  
  // END OF THE CODE.
  PrepareEnd(StartAll);
  return 0;
}
