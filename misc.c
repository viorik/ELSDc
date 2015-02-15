/*------------------------------------------------------------------------------

  Copyright (c) 2012-2014 viorica patraucean (vpatrauc@gmail.com)
  
  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as
  published by the Free Software Foundation, either version 3 of the
  License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU Affero General Public License for more details.

  You should have received a copy of the GNU Affero General Public License
  along with this program. If not, see <http://www.gnu.org/licenses/>.


  misc.c - This file belongs to ELSDc project (Ellipse and Line Segment 
           Detector with continuous validation).
         - It contains the definition of various functions needed in different 
           parts of the project.

------------------------------------------------------------------------------*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include "misc.h"


/*----------------------------------------------------------------------------*/
/** Error or exception handling: print error message and exit.
 */
void error( char *msg )
{
  fprintf( stderr,"%s\n", msg );
  exit(EXIT_FAILURE);
}


/*----------------------------------------------------------------------------*/
/** Check two double numbers for equality. 
 */
#define RELATIVE_ERROR_FACTOR 100.0
int double_equal( double a, double b )
{
  double abs_diff, aa, bb, abs_max;
  if( a == b ) return TRUE;

  abs_diff = fabs( a-b );
  aa = fabs(a);
  bb = fabs(b);
  abs_max = aa > bb ? aa : bb;

  /* DBL_MIN is the smallest normalized number, thus, the smallest
   number whose relative error is bounded by DBL_EPSILON. For
   smaller numbers, the same quantization steps as for DBL_MIN
   are used. Then, for smaller numbers, a meaningful "relative"
   error should be computed by dividing the difference by DBL_MIN. */
  if( abs_max < DBL_MIN ) abs_max = DBL_MIN;

  /* equal if relative error <= factor x eps */
  return ( abs_diff / abs_max ) <= (RELATIVE_ERROR_FACTOR * DBL_EPSILON);
}


/*----------------------------------------------------------------------------*/
/** Return +1 for positive value, and -1 for strictly negative value.
 */
int sign( double val )
{
  if( val >= 0 ) return 1;
  else return -1;
}


/*----------------------------------------------------------------------------*/
/** Euclidean distance between two points (x1,y1) and (x2,y2).
 */
double dist( double x1, double y1, double x2, double y2 )
{
  return sqrt( (x2-x1) * (x2-x1) + (y2-y1) * (y2-y1) );
}


/*----------------------------------------------------------------------------*/
/** Absolute value angle difference.
 */
double angle_diff( double a, double b )
{
  a -= b;
  while( a <= -M_PI ) a += M_2__PI;
  while( a >   M_PI ) a -= M_2__PI;
  if( a < 0.0 ) a = -a;
  return a;
}


/*----------------------------------------------------------------------------*/
/** Signed angle difference.
 */
double angle_diff_signed( double a, double b )
{
  a -= b;
  while( a <= -M_PI ) a += M_2__PI;
  while( a >   M_PI ) a -= M_2__PI;
  return a;
}


/*----------------------------------------------------------------------------*/
/** Swap two integer values.
 */
void swap_int( int *a, int *b )
{
  int t;
  t = *a; *a = *b; *b = t;
}


/*----------------------------------------------------------------------------*/
/** Swap two double values.
 */
void swap_double( double *a, double *b )
{
  double t;
  t = *a; *a = *b; *b = t;
}


/*----------------------------------------------------------------------------*/
/** Swap pointers' addresses for two arrays of points.
 */
void swap_ptr_pts( Point **a, Point **b )
{
  Point *p;

  /* check parameters */
  if( *a == NULL ) error("swap_ptr_pts: invalid input array.");
  if( *b == NULL ) error("swap_ptr_pts: invalid input array.");
   p = *a;
  *a = *b;
  *b =  p;
}


/*----------------------------------------------------------------------------*/
/** Compute the min value of an array of doubles.
 */
double min_array( double *v, int sz )
{
  double m = DBL_MAX;
  int i;

  /* check parameters */
  if( v == NULL ) error("min_array: null input array.");
  if( sz <= 0 ) error("min_array: size must be strictly positive.");

  for( i=0; i<sz; i++ ) 
    if( v[i] < m ) 
      m = v[i];

  return m;
}


/*----------------------------------------------------------------------------*/
/** Compute the min value and its index in an array of doubles. 
 */
double min_array_pos( double *v, int sz, int *pos )
{
  double m = DBL_MAX;
  int i;

  /* check parameters */
  if( v == NULL ) error("min_array_pos: null input array.");
  if( sz <= 0 ) error("min_array_pos: size must be strictly positive.");
  if( pos == NULL) error("min_array_pos: pos must be non-null.");

  for( i=0; i<sz; i++ ) 
    if( v[i] < m ) 
      {
        m = v[i];
        *pos = i;
      }
  return m;
}


/*----------------------------------------------------------------------------*/
/** Compute max element in an array and return the max value and its position.
 */
double max_array_pos( double *a, int sz, int *poz )
{
  int i;
  double max = (double)LONG_MIN;

  /* check parameters */
  if (a == NULL) error("max_array_pos: invalid input array.");
  
  for( i=0; i<sz; i++ )
    if( max<a[i] ) 
      {
        max = a[i]; 
        *poz = i;
      }
  return max;
}


/*----------------------------------------------------------------------------*/
/** Check if a value is in a given range. Return TRUE or FALSE. 
 */
int in_interval( double x, double a, double b )
{
  if( x>=a && x<b ) return TRUE;
  return FALSE;
}


/*----------------------------------------------------------------------------*/
/** Are two angles similar up to precision 'prec'?
 */
int is_aligned( double a, double theta, double prec )
{
  /* check parameters */
  if( prec < 0.0 ) error("is_similar: 'prec' must be positive.");

  /* if angle not defined, it is considered NON-similar */
  if( a == NOTDEF ) return FALSE;  /* there is no need to call the function
                                      'double_equal' here because there is
                                      no risk of problems related to the
                                      comparison doubles, we are only
                                      interested in the exact NOTDEF value */

  /* it is assumed that 'theta' and 'a' are in the range [-pi,pi] */
  theta -= a;
  if( theta < 0.0 ) theta = -theta;
  if( theta > M_3_2_PI )
    {
      theta -= M_2__PI;
      if( theta < 0.0 ) theta = -theta;
    }

  return theta <= prec;
}


/*----------------------------------------------------------------------------*/
/** Normalised difference between two angles. 
 */
double norm_angle_diff( double a, double b )
{
  /* check if valid angles */
  if( (a == NOTDEF) || (b == NOTDEF) ) return 1.0;

  /* normalize a */
  while( a <= -M_PI ) a += M_2__PI;
  while( a >   M_PI ) a -= M_2__PI;

  /* normalize b */
  while( b <= -M_PI ) b += M_2__PI;
  while( b >   M_PI ) b -= M_2__PI;

  /* compute difference and return normalised value */
  a -= b;
  while( a <= -M_PI ) a += M_2__PI;
  while( a >   M_PI ) a -= M_2__PI;

  return fabs(a) / M_PI;
}
