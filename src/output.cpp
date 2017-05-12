#include <healpix_map.h>
#include <healpix_map_fitsio.h>
#include <alm.h>
#include <fstream>
#include <iomanip>
#include "definitions.h"
#include "ParameterList.hpp"
#include "Utilities.hpp"


// Prints one single map to FITS and/or TGA file.
void WriteMap(const Healpix_Map<MAP_PRECISION> & map, std::string filename, bool inform) {
  if (filename!="0") {
    // Write to FITS:
    write_Healpix_map_to_fits("!"+filename, map, planckType<MAP_PRECISION>()); // Filename prefixed by ! to overwrite.
    if(inform==1) std::cout << ">> Healpix map written to "<<filename<<std::endl;
  } 
}


// Prints one single alm table to a TEXT file.
void WriteAlm(const Alm<xcomplex <ALM_PRECISION> > & a, const ParameterList & config, std::string keyword, bool inform) {
  std::string filename;
  std::ofstream outfile; 
  int lminout, lmaxout, mmax, l, m;

  // If requested, write alm's to file:
  if (config.reads(keyword)!="0") {
    filename = config.reads(keyword);
    outfile.open(filename.c_str());
    if (!outfile.is_open()) warning("WriteAlm: cannot open "+filename+" file.");
    outfile << "# l, m, Re(alm), Im(alm)"<<std::endl<<std::endl;
    
    // Check if output bounds are correct:
    lminout = 0;
    lmaxout = config.readi("LMAX");
    if (lmaxout > a.Lmax()) { 
      lmaxout = a.Lmax(); 
      warning("WriteAlm: LMAX beyond available data, will use the latter instead.");
    }
    mmax = -1;
    if (mmax>lminout) error ("WriteAlm: current code only allows MMAX_OUT <= LRANGE_OUT lower bound.");

    // Output all alm's:
    if (mmax<0) {
      for(l=lminout; l<=lmaxout; l++)
	for(m=0; m<=l; m++) {
	  outfile << l <<" "<< m;
#if USEXCOMPLEX
	  outfile <<" "<<std::setprecision(10)<< a(l,m).re<<" "<<std::setprecision(10)<< a(l,m).im;
#else
	  outfile <<" "<<std::setprecision(10)<< a(l,m).real()<<" "<<std::setprecision(10)<< a(l,m).imag();
#endif
	  outfile<<std::endl;
	} 
    }
    // Truncate m in alm output:
    else {
     for(l=lminout; l<=lmaxout; l++)
	for(m=0; m<=mmax; m++) {
	  outfile << l <<" "<< m;
#if USEXCOMPLEX
	  outfile <<" "<<std::setprecision(10)<< a(l,m).re<<" "<<std::setprecision(10)<< a(l,m).im;
#else
	  outfile <<" "<<std::setprecision(10)<< a(l,m).real()<<" "<<std::setprecision(10)<< a(l,m).imag();
#endif
	  outfile<<std::endl;
	}  
    }
    outfile.close();
    if(inform==1) std::cout << ">> "+keyword+" written to "+filename<<std::endl;
  }
}
