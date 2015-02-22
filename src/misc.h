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

  misc.h - This file belongs to ELSDc project (Ellipse and Line Segment 
           Detector with continuous validation).
         - It defines constants, data types, prototypes for misc functions 
           needed in different parts of the project. 

------------------------------------------------------------------------------*/


#ifndef MISC_H
#define MISC_H

/*----------------------------------------------------------------------------*/
/** Useful constants and macros.
 */

#ifndef M_LN10
#define M_LN10 2.30258509299404568402
#endif /* !M_LN10 */

#ifndef M_PI
#define M_PI   3.14159265358979323846
#endif /* !M_PI */

#ifndef FALSE
#define FALSE 0
#endif /* !FALSE */

#ifndef TRUE
#define TRUE 1
#endif /* !TRUE */

#define NOTDEF         -1024.0
#define M_3_2_PI       4.71238898038
#define M_1_2_PI       1.57079632679
#define M_2__PI        6.28318530718
#define SQRT2          1.41421356237
#define LOG10_2        0.30102999566
#define LOG10_9        0.95424250944

#define NOTUSED 0
#define USED 1

#define max(A, B) (((A)>(B))?(A):(B))
#define min(A, B) (((A)<(B))?(A):(B))


/*----------------------------------------------------------------------------*/
/** Integer cartesian coordinates of a point.
 */
typedef struct 
{
  int x, y;
} Point;

/*----------------------------------------------------------------------------*/
/** Double cartesian coordinates of a point.
 */
typedef struct 
{
  double x, y;
} PointD;


/* various function prototypes */
void error( char *msg );
int double_equal( double a, double b);
int sign( double val );
double dist( double x1, double y1, double x2, double y2 );
double angle_diff( double a, double b );
double angle_diff_signed( double a, double b );
void swap_int( int *a, int *b );
void swap_double( double *a, double *b );
void swap_ptr_pts( Point **a, Point **b );
double min_array( double *v, int sz );
double min_array_pos( double *v, int sz, int *pos );
double max_array_pos( double *a, int sz, int *poz );
int in_interval( double x, double a, double b );
int is_aligned( double a, double theta, double prec );
double norm_angle_diff( double a, double b );

#endif
