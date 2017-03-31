#include <healpix_map.h>
#include "definitions.h"

// Goes from completeness map to density contrast noise map:
void Completeness2NoiseMap(Healpix_Map<MAP_PRECISION> & Map);

// Creates binary mask from non-binary map:
void NoiseMap2BinaryMask(Healpix_Map<MAP_PRECISION> & Map);
