#include "definitions.h"
#include "ParameterList.hpp"

// Prints one single map to FITS and/or TGA file.
void WriteMap(const Healpix_Map<MAP_PRECISION> & map, std::string filename, bool inform) {
  std::string filename, tgafile;
  char *arg[4];
  char message1[100], message2[100];
  char opt1[]="-bar";

  if (filename!="0") {
    // Write to FITS:
    filename=config.reads(keyword); 
    write_Healpix_map_to_fits("!"+filename, map, planckType<MAP_PRECISION>()); // Filename prefixed by ! to overwrite.
    if(inform==1) std::cout << ">> "<<keyword<<" written to "<<filename<<std::endl;
  } 
}
