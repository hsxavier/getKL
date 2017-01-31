#include <cstring>
#include <fitsio.h>
#include "Utilities.hpp"
#include <arr.h>         // For PrepRingWeights function.


//#define HEALPIX_DATA is defined during compilation.

// Reads a Healpix FITS file containing the map2alm weights into a double array:
int ReadHealpixData(int nside, int col, double *weights, int datatype) {
  char message[200];
  std::string filename;
  fitsfile *fpointer;
  int status=0, anynul=0;
  long i, firstrow=1, firstelem=1, nelements;
  double *nulval;
  
  filename.assign(HEALPIX_DATA);
  if (filename.at(filename.length()-1)!='/') filename = filename+"/";
  // Find out what kind of data will be loaded:
  switch (datatype) {
  case 1:
    // Load map2alm weights:
    filename = filename+"weight_ring_n"+ZeroPad(nside, 10000)+".fits";
    nelements=2*nside;
    break;
  case 2:
    // Load pixel window functions:
    filename = filename+"pixel_window_n"+ZeroPad(nside, 1000)+".fits";
    nelements=4*nside+1;
    break;
  default:
    error("ReadHealpixData: unknown Healpix data type.");
  }

  // Open file:
  fits_open_table(&fpointer, filename.c_str(), READONLY, &status);
  if (status!=0) {
    sprintf(message, "ReadHealpixData: could not open table in FITS file, ERR=%d", status);
    warning(message);
  }
  
  // Prepare to, read and check for errors:
  nulval = vector<double>(0, nelements-1);
  for(i=0; i<nelements; i++) nulval[i]=666.0;
  fits_read_col(fpointer, TDOUBLE, col, firstrow, firstelem, nelements, nulval, weights, &anynul, &status);
  if (status!=0) {
    sprintf(message, "ReadHealpixData: problem reading column in FITS file table, ERR=%d", status);
    warning(message);
  }
  if(anynul!=0) {
    warning("ReadHealpixData: found NULL values in FITS table");
    printf("They are:\n");
    for (i=0; i<nelements; i++) printf("%g ",nulval[i]);
    printf("\n");
  }
  free_vector(nulval, 0, nelements-1);

  // Close file and exit:
  fits_close_file(fpointer, &status);
  if (status!=0) {
    sprintf(message, "ReadHealpixData: could not close FITS file, ERR=%d", status);
    warning(message);
  }
  return status;
}


// Prepare weights used by map2alm. weight array must be allocated already:
void PrepRingWeights(int col, arr<double> & weight, int nside) {
  double *tempweight;
  int status;
  long i;

  Announce("   Loading Healpix map weights... ");
  tempweight = vector<double>(0, 2*nside-1);
  status = ReadHealpixData(nside, col, tempweight, 1);
  if (status==0) for (i=0; i<2*nside; i++) weight[i]=1.0+tempweight[i];
  else { 
    warning("PrepRingWeights: could not load Healpix weights, using 1.0 instead.");
    weight.fill(1);
  }
  free_vector(tempweight, 0, 2*nside-1);
  Announce();
}
