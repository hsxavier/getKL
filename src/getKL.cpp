/* Code for computing the Karhunen-Loeve (KL) modes for a galaxy survey 
   with separable selection function. The idea is that we first compute 
   spherical logarithmic waves modes (angular part are expanded in 
   spherical harmonics and the radial part are logarithmic waves), and 
   then we find the eigenvalues and eigenvectors of the covariance matrix 
   to build a reduced set of linear combinations of these modes that best 
   describe the data, thus obtaining the KL modes. 
*/

#include <iostream>
#include <iomanip>               // For std::setprecision.
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
#include "definitions.h"
#include "HealpixMapProcess.hpp"
#include "Cosmology.hpp"
#include "radialTrafos.hpp"
#include "output.hpp"



/*********************************************************************************/
/***                                 Main code                                 ***/
/*** KL modes from theoretical P(k) and radial and angular selection functions ***/
/*********************************************************************************/

int main(int argc, char *argv[]) {
  using std::cout;
  using std::endl;
  using namespace ParDef; ParameterList config;         // Easy configuration file use.
  // Special variables:
  Cosmology cosmo;                                      // Cosmological parameters.
  time_t StartAll;                                      // For timing the code run.
  Healpix_Map<MAP_PRECISION> NoiseMap;
  Alm<xcomplex <ALM_PRECISION> > Nlm, Wlm;
  CovMatrix noise, angular, constant;
  int Nside, Scheme, qmax, lmax, CovN, AngN, Nfft;
  // Generic variables and iterators:
  std::string ExitAt, str1;
  std::ofstream outfile;
  int i, j, k, q, Q, l, m, L, M, ll, mm;
  double *wrapper[2], SelecScale;
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
  cosmo.load(config);


  /********************************************************/
  /*** Part 1: Load and process selection function data ***/
  /********************************************************/


  /*** Part 1.1: Prepare radial data ***/
  cout << "** Prepare radial data:"<<endl;

  // Load redshift distribution:
  Announce("   Loading redshift distribution...");
  double *rDist, *GalDens;
  LoadVecs(wrapper, config.reads("Z_DIST"), &k, &i);
  if (i!=2) warning("gelKL: Expecting two columns, found a different number.");
  rDist   = wrapper[0];
  GalDens = wrapper[1];
  Announce();

  // Re-scale redshift distribution according to SELEC_SCALE:
  Announce("   Re-scaling z dist. by SELEC_SCALE...");
  SelecScale = config.readd("SELEC_SCALE");
  for (i=0; i<k; i++) GalDens[i] = SelecScale * GalDens[i];
  Announce();

  // Compute density, output result and/or exit if requested: 
  Announce("   Converting to galaxy density in gals/(h^-1 Mpc)^3...");
  SphDens2CartDens(cosmo, rDist, GalDens, k);
  Announce();
  OutputTable(config, "GALDENS_OUT", wrapper, k, 2);
  if (ExitAt=="GALDENS_OUT") { PrepareEnd(StartAll); return 0; }
  
  // Compute noise term, output result and/or exit if requested:
  Announce("   Compute radial noise term 1/n(r)...");
  for(i=0; i<k; i++) {
    if (GalDens[i]<0) warning("getKL: Found negative galaxy density.");
    if (GalDens[i]>0) GalDens[i] = 1.0/GalDens[i];
  }
  Announce();
  OutputTable(config, "RADNOISE_OUT", wrapper, k, 2);
  if (ExitAt=="RADNOISE_OUT") { PrepareEnd(StartAll); return 0; }

  // Compute logarithmic wave transform of radial part of the noise:
  std::vector< std::complex <double> > ZZr2trafo;
  double rStart, rEnd, rWin0, rWin1;
  qmax   = config.readi("QMAX"); // This is max. q for the radial modes.
  Nfft   = 4*qmax;
  rStart = ComDist(cosmo, config.readd("TRAFO_ZRANGE",0));
  rEnd   = ComDist(cosmo, config.readd("TRAFO_ZRANGE",1));
  rWin0  = ComDist(cosmo, config.readd(  "SEL_ZRANGE",0));
  rWin1  = ComDist(cosmo, config.readd(  "SEL_ZRANGE",1));

  cout << "   Trafo from "<<rStart<<" to "<<rEnd<<", with non-zero values at "<<rWin0<< " < r < "<<rWin1<<" h^-1 Mpc."<<endl;
  cout << "   Will compute covariance for "<<2*qmax+1<<" radial modes."<<endl;
  i = config.readi("ZZR2_METHOD");
  if (i!=1 && i!=2) { warning("Unknown ZZR2_METHOD, will use method 1."); i=1; }
  // Use FFT to compute radial noise trafo:
  if (i==1) {
    Announce("   Computing trafo with FFTW...");
    ZZr2TranformFFT(rDist, GalDens, k, rWin0, rWin1, rStart, rEnd, Nfft, ZZr2trafo);
  }
  // Use Romberg Integration to compute radial noise trafo:
  if (i==2) {
    Announce("   Computing trafo with Romberg...");
    ZZr2TranformRomb(rDist, GalDens, k, rWin0, rWin1, rStart, rEnd, Nfft, ZZr2trafo);
  }
  Announce();
  // Output result if requested:
  if (config.reads("ZZR2_OUT")!="0") {
    str1 = config.reads("ZZR2_OUT");
    Announce(">> Writing ZZR2_OUT to "+str1+"...");
    PrintZZr2(str1, rStart, rEnd, Nfft, ZZr2trafo);
    Announce();
  }
  if (ExitAt=="ZZR2_OUT") { PrepareEnd(StartAll); return 0; }


  /*** Part 1.2: Prepare angular data ***/
  cout << "** Prepare angular data:"<<endl;

  // Read in completeness map and compute density contrast noise map:
  Announce("   Loading completeness map...");
  read_Healpix_map_from_fits(config.reads("COMPLETE_MAP"), NoiseMap);
  Nside  = NoiseMap.Nside();
  Scheme = NoiseMap.Scheme();
  if (Scheme !=RING) warning("Expecting map in the RING ordering scheme.");
  Announce();
  cout << "   COMPLETE_MAP has Nside = "<<Nside<<".\n";
  Announce("   Computing density contrast noise map...");
  Completeness2NoiseMap(NoiseMap);
  Announce();
  // Output noise map if requested:
  if (config.reads("NOISEMAP_OUT")!="0") {
    str1 = config.reads("NOISEMAP_OUT");
    Announce(">> Writing NOISEMAP_OUT to "+str1);
    write_Healpix_map_to_fits("!"+str1, NoiseMap, planckType<MAP_PRECISION>());
    Announce();
  }
  if (ExitAt=="NOISEMAP_OUT") { PrepareEnd(StartAll); return 0; }

  // Get harmonic coefficients from noise map:
  arr<double> RingWeights(2*Nside);
  PrepRingWeights(1, RingWeights, Nside);
  Announce("   Get harmonic coefficients from noise map...");
  Nlm.Set(2*lmax,2*lmax); // Noise cov. matrix depends on integral of 3 Ylm up to Nlm with l=2*LMAX.
  for(l=0; l<=2*lmax; l++) for (m=0; m<=l; m++) { Nlm(l,m).real(0); Nlm(l,m).imag(0); }
  map2alm(NoiseMap, Nlm, RingWeights);
  Announce();


  /*** Part 1.3: Compute i_wlm i*_LWM, the constant term of cov. matrix if the mean is not subtracted ***/
  
  // Create binary angular mask:
  Announce("   Creating binary angular mask...");
  NoiseMap2BinaryMask(NoiseMap);
  Announce();
  // Output noise map if requested:
  if (config.reads("MASK_OUT")!="0") {
    str1 = config.reads("MASK_OUT");
    Announce(">> Writing MASK_OUT to "+str1);
    write_Healpix_map_to_fits("!"+str1, NoiseMap, planckType<MAP_PRECISION>());
    Announce();
  }
  if (ExitAt=="MASK_OUT") { PrepareEnd(StartAll); return 0; }
  
  // Compute harmonic coefficients of the mask: 
  Announce("   Get harmonic coefficients from binary mask...");
  Wlm.Set(lmax,lmax);
  for(l=0; l<=lmax; l++) for (m=0; m<=l; m++) { Wlm(l,m).real(0); Wlm(l,m).imag(0); }
  map2alm(NoiseMap, Wlm, RingWeights);
  // Deallocate memory that will not be used anymore:
  NoiseMap.SetNside(1, RING);
  Announce();
  WriteAlm(Wlm, config, "MASK_ALM_OUT", 1);
  if (ExitAt=="MASK_ALM_OUT") { PrepareEnd(StartAll); return 0; }

  /*
  // Compute radial transform:
  std::vector< std::complex<double> > intZr2;
  double U;
  U = IntervalU(rStart, rEnd);
  intZr2.resize(2*qmax+1);
  for (q=-qmax; q<=qmax; q++) {
    z.real(1.5);  z.imag(q2w(q,U));
    intZr2[q+qmax] = (ZetawConj(q2w(q,U),rWin1) - ZetawConj(q2w(q,U),rWin0)) / z;  // Compared with NIntegrate in mathematica, good.
  }
  
  // Create constant part of cov. matrix:
  constant.Alloc(qmax,lmax);
  for (q=-qmax; q<=qmax; q++) for(l=0; l<=lmax; l++) for (m=-l; m<=l; m++) {
	for (Q=-qmax; Q<=qmax; Q++) for(L=0; L<=lmax; L++) for (M=-L; M<=L; M++) {
	      if (m<0) z = MinusOneToPower(m)*conj(Wlm(l,-m));
	      else z = Wlm(l,m);
	      if (M<0) z2 = MinusOneToPower(M)*conj(Wlm(L,-M));
	      else z2 = Wlm(L,M);	      
	      z = intZr2[q+qmax]*z * conj(intZr2[Q+qmax]*z2);
	      //z = intZr2[q+qmax] * conj(intZr2[Q+qmax]);
	      //z = z * z2;
	      constant.set(q,l,m,Q,L,M, z);
	    }
      }
  
  // If requested, print out the covariance matrix:
  str1=config.reads("COVCONST_OUT");
  if (str1!="0") {
    Announce(">> Writting constant part of Cov. Matrix to "+str1+"...");
    constant.Print(str1, CovOut);
    Announce();
  }
  // Exit if this is the last output requested:
  if (ExitAt=="COVCONST_OUT") {
    PrepareEnd(StartAll); return 0;
  }
  */


  /***********************************************/
  /*** Part 2: Compute noise covariance matrix ***/
  /***********************************************/


  /*** Part 2.1: Compute angular part of the noise covariance matrix ***/
  cout << "** Prepare NOISE covariance matrix:"<<endl;

  // Prepare Covariance matrix object:
  Announce("   Allocating memory for matrix...");
  angular.Alloc(lmax);
  noise.Alloc(qmax,lmax);
  CovN = noise.Nentries();
  AngN = angular.Nentries();
  Announce();
  // DEBUG: print index asignment for covariance matrix:
  //printf("%6s %6s %6s %6s %6s\n", "i", "w", "l", "m", "i");
  //for (i=0; i<CovN; i++) {
  //  noise.i2qlm(i,&k,&l,&m);
  //  printf("%6d %6d %6d %6d %6d\n", i, k, l, m, noise.qlm2i(k,l,m));
  //}
  //return 0;

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
    // Compute the element given by a sum over Nlm:
    z.real(0); z.imag(0);
    for (ll=FirstEll(L,l,M,m); ll<=L+l; ll=ll+2) {
      if (M>m) z2 = MinusOneToPower(M-m) * conj(Nlm(ll,M-m));
      else     z2 = Nlm(ll,m-M);
      z += z2 * ThreeYlmIntegral(ll, m-M, L, M, l, m);
    }
    // Save element to the matrix:
    angular.set(i,j,z);
  }
  Announce();

  // If requested, print out the covariance matrix:
  str1=config.reads("ANGCOV_OUT");
  if (str1!="0") {
    Announce(">> Writting angular part of Cov. Matrix to "+str1+"...");
    angular.Print(str1, CovOut);
    Announce();
  }
  // Exit if this is the last output requested:
  if (ExitAt=="ANGCOV_OUT") {
    PrepareEnd(StartAll); return 0;
  }
 

  // Compute the noise angular covariance matrix: 
  cout <<  "   Cov. matrix side is N="<<CovN<<endl;
  Announce("   Computing diag. and below of full cov. matrix...");
  long2=(((long)CovN+1)*((long)CovN))/2;
  // LOOP over Cov. Matrix elements:
#pragma omp parallel for schedule(dynamic) private(i,j,Q,L,M,q,l,m,z,ll,z2)
  for (long1=0; long1<long2; long1++) {
    i = (int)((sqrt(8*long1+1)-1.0)/2.0);
    j = (int)(long1-(i*(i+1))/2);
    noise.i2qlm(i,&q,&l,&m);
    noise.i2qlm(j,&Q,&L,&M);
    // Multiply the angular part by the radial part:
    if (angular.qlm2i(0,l,m)>angular.qlm2i(0,L,M))
      z = angular(0,l,m,0,L,M) * ZZr2trafo[qq2i(Q-q, Nfft)];
    else 
      z = conj(angular(0,L,M,0,l,m)) * ZZr2trafo[qq2i(Q-q, Nfft)];
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
