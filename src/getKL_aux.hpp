#ifndef GETKL_AUX
#define GETKL_AUX 1

// Returns (-1)^m:
double MinusOneToPower(int m);

// Returns the integral over all sphere of: Y_l1m1 x Y_l2m2 x Y*_l3m3: 
double ThreeYlmIntegral(int l1, int m1, int l2, int m2, int l3, int m3);

#endif
