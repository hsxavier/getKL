#ifndef FITSFUNCTIONS_H
#define FITSFUNCTIONS_H 1

#include <arr.h>
 
int ReadHealpixData(int nside, int col, double *weights, int datatype);
// Prepare weights used by map2alm. weight array must be allocated already:
void PrepRingWeights(int col, arr<double> & weight, int nside);

#endif
