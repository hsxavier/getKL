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
#include <algorithm>             // for max function.
#include "getKL_aux.hpp"
#include "fitsfunctions.hpp"     // For PrepRingWeights function.
#include "Utilities.hpp"
#include "CovMatrix.hpp"
#include "ParameterList.hpp"

#define MAP_PRECISION double
#define ALM_PRECISION double



/*********************************************************************************/
/***                                 Main code                                 ***/
/*** KL modes from theoretical P(k) and radial and angular selection functions ***/
/*********************************************************************************/

int main(int argc, char *argv[]) {
  using std::cout;
  using std::endl;
  using namespace ParDef; ParameterList config;         // Easy configuration file use.
  // Special variables:
  time_t StartAll;                                      // For timing the code run.
  Healpix_Map<MAP_PRECISION> NoiseMap;
  Alm<xcomplex <ALM_PRECISION> > Nlm;
  CovMatrix noise;
  int Nside, Scheme, lmax, CovN;
  // Generic variables and iterators:
  std::string ExitAt, str1;
  int i, j, k, l, m, L, M, ll, mm;  
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
  ExitAt = config.reads("EXIT_AT");
  if (ExitAt!="0") config.findpar(ExitAt+":"); // Produce warning if last output is unknown.
  lmax = config.readi("LMAX");
  

  /********************************************************/
  /*** Part 1: Load and process selection function data ***/
  /********************************************************/
  
  // Read in noise map:
  Announce("Loading noise map...");
  read_Healpix_map_from_fits(config.reads("NOISE_MAP"), NoiseMap);
  Nside  = NoiseMap.Nside();
  Scheme = NoiseMap.Scheme();
  if (Scheme !=RING) warning("Expecting map in the RING ordering scheme.");
  Announce();
  cout << "NOISE_MAP has Nside = "<<Nside<<".\n";

  // Get harmonic coefficients from noise map:
  arr<double> RingWeights(2*Nside);
  PrepRingWeights(1, RingWeights, Nside);
  Announce("Get harmonic coefficients from noise map...");
  Nlm.Set(2*lmax,2*lmax); // Noise cov. matrix depends on integral of 3 Ylm up to Nlm with l=2*LMAX.
  for(l=0; l<=2*lmax; l++) for (m=0; m<=l; m++) { Nlm(l,m).real(0); Nlm(l,m).imag(0); }
  map2alm(NoiseMap, Nlm, RingWeights);
  Announce();


  /***********************************************/
  /*** Part 2: Compute noise covariance matrix ***/
  /***********************************************/

  /*** Part 2.1: Compute angular part of the noise covariance matrix ***/

  // Prepare Covariance matrix object:
  Announce("Allocating memory for covariance matrix...");
  noise.Alloc(lmax);
  CovN = noise.Nentries();
  Announce();
  cout << "Cov. matrix side is N="<<CovN<<endl;

  // Compute the noise angular covariance matrix: 
  Announce("Computing diag. and lower triangular of cov. matrix...");
  long2=(((long)CovN+1)*((long)CovN))/2;
  // LOOP over Cov. Matrix elements:
#pragma omp parallel for schedule(dynamic) private(i,j,k,L,M,l,m,z,ll,z2)
  for (long1=0; long1<long2; long1++) {
    i = (int)((sqrt(8*long1+1)-1.0)/2.0);
    j = (int)(long1-(i*(i+1))/2);
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


  // If requested, print out the covariance matrix:
  str1=config.reads("COVMATRIX_OUT");
  if (str1!="0") {
    Announce(">> Writting Cov. Matrix to "+str1+"...");
    noise.Print(str1, CovOut);
    Announce();
  }
  
  // Exit if this is the last output requested:
  if (ExitAt=="COVMATRIX_OUT") {
    PrepareEnd(StartAll); return 0;
  }

  


  /************************************************/
  /*** Part 3: Get eigenvalues and eigenvectors ***/
  /************************************************/

  // Find eigenvalues and eigenvectors:
  Announce("Solving eigensystem...");
  noise.EigenSolve();
  Announce();

  // Copy to external data:
  Announce("Copying solution to external data...");
  std::vector<double> Eigenvalues;
  std::vector< std::vector< std::complex<double> > > Eigenvectors;
  
  Eigenvalues.resize(CovN);
  Eigenvectors.resize(CovN);
  for(i=0; i<CovN; i++) Eigenvectors[i].resize(CovN);

  noise.Eigenvalues2(Eigenvalues);
  noise.Eigenvectors2(Eigenvectors);
  Announce();
  
  // Print solution to the screen:
  cout << endl << "Eigenvalues:\n";
  for (i=0;i<CovN; i++) cout << Eigenvalues[i] << " ";
  cout << endl;
  cout << endl << "Eigenvectors:\n";
  for (i=0; i<CovN; i++) {
    for (j=0; j<CovN; j++) cout << Eigenvectors[i][j] << " ";
    cout<<endl;
  }
  

  /**************************/
  /*** End of the program ***/
  /**************************/

  PrepareEnd(StartAll);
  return 0;
}
