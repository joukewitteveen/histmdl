#include <assert.h>
#include <limits.h>
#include <R.h>

#define likely( x ) __builtin_expect( !!( x ), 1 )


/* Highest index of a value less than or equal to a search value
 * x: sorted array of values
 * l: length of x
 * v: search value
 * i: index, initialized to 0
 */
void maxLE( double const * x, int const * l, double const * v, int * i ){
  assert( *i == 0 );
  if( *l <= 0 ) return;
  // The highest index in x
  int max = *l - 1;
  // Initialize d to the floor of the base 2 logarithm of max
  int d = max | ( max >> 1 );
  d |= d >> 2; d |= d >> 4; d |= d >> 8; d |= d >> 16;
#if INT_MAX > (1 << 32)
  d |= d >> 32;
#endif
  d ^= d >> 1;
  do {
    // Verify that *i + d is within range
    if( likely( max & d ) ){
      if( likely( x[*i | d] > *v ) ) {
        // No more range checking is needed
        while( d >>= 1 ) if( x[*i | d] <= *v ) *i |= d;
        break;
      } else *i |= d;
    } //if
  } while( d >>= 1 );
  if( *x <= *v ) ++*i;
} //maxLE

