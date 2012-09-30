#ifndef __ELOG_H__
#define __ELOG_H__

#include <math.h>

/**
 * "Extended" logarithm functions that are used to work with small
 * numbers that would normally result in float underflows. The methods
 * take into account the possibility of log(0) and return appropriate
 * values, removing the necessity to handle log(0) cases from the main
 * codebase. These functions are especially useful for the
 * multiplication of small probability values.
 *
 * These functions are implemented roughly following the pseudo code
 * from "Numerically Stable Hidden Markov Model Implementation" by
 * Tobias P. Mann.
 */

#ifndef INFINITY
#define INFINITY HUGE_VAL
#endif
#define LOG_ZERO -INFINITY


/**
 * Computes e^x, returning 0 if x == LOG_ZERO
 */
#define eexp(x) (((x) == LOG_ZERO) ? 0 : exp(x))

/**
 * Computes ln(x), returning LOG_ZERO if x==0
 */
#define eln(x) (((x) == 0) ? LOG_ZERO : log(x))

/**
 * Allows computation of ln(x+y) in a numerically stable way given ln(x)
 * and ln(y) as arguments (i.e. computes ln(e^ln(x) + e^ln(y))). By
 * exponentiating only the smaller of the two values, precision is
 * maintained when the values are much different magnitudes (and
 * overflows are avoided when possible).  Uses the following identity:
 *
 * ln(x + y) = ln(x) + ln(x + y) - ln(x)
 *           = ln(x) + ln((x +y) / x)
 *           = ln(x) + ln(1 + y/x)
 *           = ln(x) + ln(1 + e^ln(y/x))
 */
#define elnsum(lnx, lny) \
  (((lnx) == LOG_ZERO || (lny) == LOG_ZERO) ?	\
   (((lnx) == LOG_ZERO) ? (lny) : (lnx)) :		\
   (((lnx) > (lny)) ? (lnx) + eln(1 + exp((lny)-(lnx))) :	\
    (lny) + eln(1 + exp((lnx)-(lny)))))


/**
 * Allows computation of ln(x*y) in a numerically stable way, given ln(x)
 * and ln(y) as arguments (i.e. computes ln(x)+ln(y)).
 * If one of the arguments is LOG_ZERO, LOG_ZERO is returned.
 */
#define elnproduct(lnx, lny) \
  (((lnx) == LOG_ZERO) || ((lny) == LOG_ZERO)) ? LOG_ZERO : (lnx) + (lny)


#endif
