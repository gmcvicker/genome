
#include <math.h>

#include "stat.h"


/**
 * Returns value of normal cumulative distribution function evaluated at x
 */
double stat_norm_cdf(double mean, double sd, double x) {
  return 0.5 * (1.0 + erf((x-mean)/(sd*M_SQRT2)));
}
