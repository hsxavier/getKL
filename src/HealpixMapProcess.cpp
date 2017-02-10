#include <healpix_map.h>
#include "definitions.h"

/* Uses a completeness map (including the survey and veto masks and the decrease due to seeing, 
galactic extinction and star density, usually assigned to systematic weights) to compute the 
Poisson noise term for the galaxy density normalized by the mean density, n(r)/<n>. 
It transforms the input into the output. */
void Completeness2NoiseMap(Healpix_Map<MAP_PRECISION> & Map) {
  int i, Npix;
  Npix = 12*Map.Nside()*Map.Nside();
  // If masked (completeness<=0), leave it like that (noise there is zero, just like the signal);
  // If unmasked, invert it (noise there is amplified by the weighting):
  for (i=0; i<Npix; i++) if (Map[i] > 0.0) Map[i] = 1.0/Map[i]; 
}
