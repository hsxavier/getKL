#include <gsl/gsl_sf_coupling.h>
#include <cmath>

// Returns (-1)^m:
double MinusOneToPower(int m) {
  if (m%2==0) return 1.0;
  else return -1.0;
}

// Returns the integral over all sphere of: Y_l1m1 x Y_l2m2 x Y*_l3m3: 
double ThreeYlmIntegral(int l1, int m1, int l2, int m2, int l3, int m3) {
  const double FourPi = 12.56637061435917;
  int twol1, twol2, twol3;
  if ( m1+m2-m3    != 0) return 0.0;
  if ((l1+l2+l3)%2 == 1) return 0.0;
  if (l3 < fabs(l1-l2) ) return 0.0;
  if (l3 > l1+l2       ) return 0.0;
  twol1 = 2*l1;
  twol2 = 2*l2;
  twol3 = 2*l3;
  return sqrt((double)((twol1+1)*(twol2+1)*(twol3+1))/FourPi) * MinusOneToPower(m3) * 
    gsl_sf_coupling_3j(twol1, twol2, twol3, 2*m1, 2*m2, -2*m3) *
    gsl_sf_coupling_3j(twol1, twol2, twol3,    0,     0,    0);
}
