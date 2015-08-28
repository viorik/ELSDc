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


  ring.c - This file belongs to ELSDc project (Ellipse and Line Segment 
           Detector with continuous validation).
         - It contains functions to work with connected regions of pixels, 
           approximated by rings.

------------------------------------------------------------------------------*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "misc.h"
#include "ring.h"


/*----------------------------------------------------------------------------*/
/** Copy one ring structure to another.
 */
void copy_ring( Ring *in, Ring *out )
{
  /* check parameters */
  if( (in == NULL) || (out == NULL) ) 
    error("copy_ring: invalid 'in' or 'out'.");

  out->x1 = in->x1;
  out->y1 = in->y1;
  out->x2 = in->x2;
  out->y2 = in->y2;
  out->width = in->width;
  out->cx = in->cx;
  out->cy = in->cy;
  out->theta = in->theta;
  out->ax = in->ax;
  out->bx = in->bx;
  out->ang_start = in->ang_start;
  out->ang_end   = in->ang_end;
  out->wmin = in->wmin;
  out->wmax = in->wmax;
  out->full = in->full;
}


/*----------------------------------------------------------------------------*/
/** Switch the delimiting points of the ring and set the switch flag. 
 */
void swap_ring( Ring *r )
{
  /* check parameters */
  if( r == NULL ) error("swap_ring: invalid ring.");

  swap_double( &(r->x1), &(r->x2) );
  swap_double( &(r->y1), &(r->y2) );
  swap_double( &(r->ang_start), &(r->ang_end) );
  r->full = 1;
}

/*----------------------------------------------------------------------------*/
/** Check if circle is valid: radius must be positive.
 */
int check_circ_ring( Ring *cring )
{
  /* check parameters */
  if( cring == NULL ) error("check_circ_ring: invalid ring.");

  /* reject if degenerate ellipse */
  if( cring->wmin*cring->wmax > 0 ) return FALSE;
  return TRUE;
}


/*----------------------------------------------------------------------------*/
/** Check if ellipse ring is valid: wmin and wmax have different signs and 
    have reasonably small values.
 */
int check_ell_ring( Ring *ering )
{
  /* check parameters */
  if( ering == NULL ) error("check_ell_ring: invalid ring.");

  /* reject if degenerate ellipse */
  if( ering->wmin*ering->wmax > 0 ) return FALSE;
  return TRUE;
}



/*----------------------------------------------------------------------------*/
/** Check if circle is valid: radius must be positive.
 */
int check_circ( double *cparam )
{
  /* check parameters */
  if( cparam == NULL ) error("check_circ: invalid circle.");

  /* reject if degenerate circle (negative or close to 0 radius) */
  if( cparam[2] <= 1 ) return FALSE;
  /* reject directly if very small arc, to win time */
  /* if( angle_diff( cring->ang_end, cring->ang_start) < 0.14 ) return FALSE;*/
  return TRUE;
}


/*----------------------------------------------------------------------------*/
/** Check if ellipse is valid: axes must be positive; ensure big axis is first 
    and orientation in [-pi/2, pi/2].
 */
int check_ell( double *eparam )
{
  /* check parameters */
  if( eparam == NULL ) error("check_ell: invalid ellipse.");

  /* reject if degenerate ellipse (negative axes) */
  if( (eparam[2] <= 1) || (eparam[3] <= 1) ) return FALSE;
  /* make sure 'ax' is the greater axis */
  if( eparam[2] < eparam[3] )
    {
      swap_double( &(eparam[2]), &(eparam[3]) );
      eparam[4] += M_1_2_PI;
    }
  /* reject directly if very thin ellipse, to win time */
  if( eparam[2] / eparam[3] > 100 ) return FALSE;

  /* make sure orientation belongs to [-pi/2, pi/2] */
  if( eparam[4] > M_1_2_PI ) eparam[4] -= M_PI;
  if( eparam[4] > M_PI ) eparam[4] -= M_PI;
  if( eparam[4] < -M_1_2_PI ) eparam[4] += M_PI;

  return TRUE;
}


/*----------------------------------------------------------------------------*/
/** Test if a given point is inside a circular ring: 
    (x,y) must be at a distance between [dmin,dmax] wrt centre at an angle 
    between ['ang_start','ang_end']. All angles are in [0 , 2pi].
 */
int is_in_circ_ring( Ring * r, int x, int y )
{
  /* check parameters */
  if( r == NULL ) error("is_in_ring: invalid input ring.");

  double d = dist( x, y, r->cx, r->cy) - r->ax;

  double theta = atan2( (double)y - r->cy, (double)x - r->cx);

  /* make sure angle in [0, 2pi] */
  if( theta<0 ) theta += M_2__PI;

  if( r->ang_end > r->ang_start ) 
    {
      if( (theta < r->ang_start) || (theta > r->ang_end) ) return FALSE;
    }
  else
    if( (theta < r->ang_start) && (theta > r->ang_end) ) return FALSE;

  if( (d < r->wmin) || (d > r->wmax) ) return FALSE;
  
  return TRUE;
}


/*----------------------------------------------------------------------------*/
/** Test if a given point is inside an elliptical ring: 
    (x,y) must be at a distance between [dmin,dmax] wrt centre at an angle 
    between ['ang_start','ang_end']. All angles are expressed in [0 , 2pi].
 */
int is_in_ell_ring( Ring * r, int x, int y )
{
  /* check parameters */
  if( r == NULL ) error("is_in_ring: invalid input ring.");

  double d = d_rosin( r, (double)x, (double)y );

  double theta = atan2( (double)y - r->cy, (double)x - r->cx);

  /* make sure angle in [0, 2pi] */
  if( theta < 0 ) theta += M_2__PI;
  
  if( r->ang_end > r->ang_start ) 
    {
      if( (theta < r->ang_start) || (theta > r->ang_end) ) return FALSE;
    }
  else
    if( (theta < r->ang_start) && (theta > r->ang_end) ) return FALSE;

  if( (d < r->wmin) || (d > r->wmax) ) return FALSE;
  return TRUE;
}


/*----------------------------------------------------------------------------*/
/** Compute the point belonging to an ellipse, closest to a given point with
    integer coordinates.
 */
void rosin_point( Ring *ering, int x, int y, double *xi, double *yi )
{  
  double xtmp, ytmp;
  double ae2, be2, fe2; 
  double xp, yp, xp2, yp2, xx, yy;  
  double delta, A, ah, bh2, term;
  double d[4];
  int pos;
  xtmp = x - ering->cx;
  ytmp = y - ering->cy;
  ae2  = ering->ax * ering->ax;
  be2  = ering->bx * ering->bx;
  fe2 = ae2 - be2;
  xp = xtmp * cos(-ering->theta) - ytmp * sin(-ering->theta);
  yp = xtmp * sin(-ering->theta) + ytmp * cos(-ering->theta);
  xp2 = xp * xp;
  yp2 = yp * yp;
  delta = (xp2 + yp2 + fe2) * (xp2 + yp2 + fe2) - 4 * xp2 * fe2;
  A = (xp2 + yp2 + fe2 - sqrt(delta))/2.0; 
  ah = sqrt(A);
  bh2 = fe2 - A;
  term = (A * be2 + ae2 * bh2);
  xx = ah * sqrt( ae2 * (be2 + bh2) / term );
  yy = ering->bx * sqrt( bh2 * (ae2 - A) / term );

  d[0] = dist( xp, yp, xx, yy );
  d[1] = dist( xp, yp, xx,-yy );
  d[2] = dist( xp, yp,-xx, yy );
  d[3] = dist( xp, yp,-xx,-yy );

  min_array_pos( d, 4, &pos );
  switch(pos)
    { 
      case 0: break;
      case 1: yy = -yy; break;
      case 2: xx = -xx; break;
      case 3: xx = -xx; yy = -yy; break;
      default: break;
    }
  *xi = xx * cos(ering->theta) - yy * sin(ering->theta) + ering->cx;
  *yi = xx * sin(ering->theta) + yy * cos(ering->theta) + ering->cy;
}


/*----------------------------------------------------------------------------*/
/** Approximate the distance between a point and an ellipse using Rosin 
    distance.
 */
double d_rosin( Ring *ering, double x, double y )
{ 
  double xtmp, ytmp;
  double ae2, be2, fe2; 
  double xp, yp, xp2, yp2, xx, yy;  
  double delta, A, ah, bh2, term;
  double d[4];

  xtmp = x - ering->cx;
  ytmp = y - ering->cy;
  ae2  = ering->ax * ering->ax;
  be2  = ering->bx * ering->bx;
  fe2 = ae2 - be2;
  xp = xtmp * cos(-ering->theta) - ytmp * sin(-ering->theta);
  yp = xtmp * sin(-ering->theta) + ytmp * cos(-ering->theta);
  xp2 = xp * xp;
  yp2 = yp * yp;
  delta = (xp2 + yp2 + fe2) * (xp2 + yp2 + fe2) - 4 * xp2 * fe2;
  A = (xp2 + yp2 + fe2 - sqrt(delta))/2.0; 
  ah = sqrt(A);
  bh2 = fe2 - A;
  term = (A * be2 + ae2 * bh2);
  xx = ah * sqrt( ae2 * (be2 + bh2) / term );
  yy = ering->bx * sqrt( bh2 * (ae2 - A) / term );

  d[0] = dist( xp, yp, xx, yy);
  d[1] = dist( xp, yp, xx,-yy);
  d[2] = dist( xp, yp,-xx, yy);
  d[3] = dist( xp, yp,-xx,-yy);
  double dmin = min_array(d,4); 
  if( xp2+yp2 > xx*xx+yy*yy ) return dmin;
  else return -dmin; 
}


/*----------------------------------------------------------------------------*/
/** Compute ellipse foci, given ellipse params.
 */
void ellipse_foci( double *param, double *foci )
{
  double f;
  /* check parameters */
  if( param == NULL) error("ellipse_foci: invalid input ellipse.");
  if( foci == NULL) error("ellipse_foci: 'foci' must be non null.");

  f = sqrt( param[2] * param[2] - param[3] * param[3] );
  foci[0] = param[0] + f * cos(param[4]);
  foci[1] = param[1] + f * sin(param[4]);
  foci[2] = param[0] - f * cos(param[4]);
  foci[3] = param[1] - f * sin(param[4]);
}


/*----------------------------------------------------------------------------*/
/** Compute the angle of the normal to a point belonging to an ellipse 
    using the focal property.
 */
double ellipse_normal_angle( double x, double y, double *foci )
{
  double tmp1, tmp2, tmp3, theta;
  /* check parameters */
  if( foci == NULL ) error("ellipse_normal_angle: 'foci' must be non null");

  tmp1 = atan2( y-foci[1], x-foci[0]);
  tmp2 = atan2( y-foci[3], x-foci[2]);
  tmp3 = angle_diff_signed( tmp1, tmp2 );

  theta = tmp1 - tmp3/2.0;
  while( theta <= -M_PI ) theta += M_2__PI;
  while( theta >   M_PI ) theta -= M_2__PI;
  return theta;
}


/*----------------------------------------------------------------------------*/
/** Determine if the gradient converges or diverges to/from the centre.
    Return 1 if gradient diverges, 0 if it converges.
    Circle case.
 */
static int dir_gradient_circ( int px, int py, double ang, double *param )
{
  double a;
  int dir = 0;
  /* check parameters */
  if( param == NULL ) error("dir_gradient_circ: invalid param.");

  a = atan2( py- param[1], px-param[0] );
  if( angle_diff( a, ang ) < M_1_2_PI )
    dir = 1;
  return dir;
}


/*----------------------------------------------------------------------------*/
/** Determine if the gradient converges or diverges to/from the centre.
    Return 1 if gradient diverges, 0 if it converges.
    Ellipse case.
 */
static int dir_gradient_ell( int px, int py, double ang, double *foci )
{
  double a;
  int dir = 0;
  /* check parameters */
  if( foci == NULL ) error("dir_gradient_ell: invalid ellipse.");
  a = ellipse_normal_angle( (double)px, (double)py, foci );
  if( angle_diff( a, ang ) < M_1_2_PI )
    dir = 1;
  return dir;
}


/*----------------------------------------------------------------------------*/
/** Compute the widths of a circle ring covering a set of pixels, given the
    centre of the circle and the radius.
 */
static void circ_ring_width( Point *reg, int reg_size, Ring *cring )
{
  int i;
  double w;
  double wmin = DBL_MAX;
  double wmax = -DBL_MAX;
  /* check parameters */
  if( reg == NULL ) error("circ_ring_width: invalid region.");
  if( cring == NULL ) error("circ_ring_width: invalid ring.");
  
  /* compute widths */ 
  for( i=0; i<reg_size; i++ )
    {  
      w = dist( reg[i].x, reg[i].y, cring->cx, cring->cy) - cring->ax;    
      if( w<wmin ) wmin = w;
      if( w>wmax ) wmax = w;
    }
  cring->wmin = wmin; cring->wmax = wmax;
}


/*----------------------------------------------------------------------------*/
/** Compute the widths of an ellipse ring covering a set of pixels, given the
    parameters of the ellipse.
 */
static void ell_ring_width( Point *reg, int reg_size, Ring *ering )
{
  int i;
  double w;
  double wmin = DBL_MAX;
  double wmax = -DBL_MAX;
  /* check parameters */
  if( reg == NULL ) error("ell_ring_width: invalid region.");
  if( ering == NULL ) error("ell_ring_width: invalid ring.");
  
  /* compute widths */ 
  for( i=0; i<reg_size; i++ )
    {  
      w = d_rosin( ering, (double)reg[i].x, (double)reg[i].y );
      if( w<wmin ) wmin = w;
      if( w>wmax ) wmax = w;
    }
  ering->wmin = wmin; ering->wmax = wmax;
}


/*----------------------------------------------------------------------------*/
/** Given conic parameters and extreme points, fill in all the ring parameters.
    If conic is 1, get circle ring, if 0 get ellipse ring.
 */
int get_ring( Point *reg, int reg_size, double ang, double *param, 
              double *pext1, double *pext2, int conic, double spir, Ring *r, 
              int *grad_dir, double *foci, int msize )
{
  double tang;
  /* check parameters */
  if( reg == NULL ) error("get_ring: invalid region."); 
  if( param == NULL ) error("get_ring: invalid conic parameters."); 
  if( r == NULL ) error("get_ring: ring must be non null.");
  if( foci == NULL ) error("get_ring: foci must be non null."); 
  
  /* fill in centre */
  r->cx = param[0]; r->cy = param[1];

  /* Init 'full' flag; by default, set to 0. */
  r->full = 0;

  if( conic == CIRCLE ) 
    {
      if( !check_circ(param) ) return FALSE;
      *grad_dir = dir_gradient_circ( reg[0].x, reg[0].y, ang, param );   
    }
  else
    {
      if( !check_ell(param) ) return FALSE;
      ellipse_foci( param, foci );
      *grad_dir = dir_gradient_ell( reg[0].x, reg[0].y, ang, foci );
    }
  /* fill in axes and orientation: this must be done after check, as for the 
     ellipse the axes might have been swapped, and orientation normalized */
  r->ax = param[2]; r->bx = param[3];
  r->theta = param[4];

  /* fill in extremal points */
  r->x1 = pext1[0]; r->y1 = pext1[1];
  r->x2 = pext2[0]; r->y2 = pext2[1]; 
  
  /* if gradient converges, interchange the extremal points to keep 
     trigonometric sense */    
  if( *grad_dir == 0 )
    {
      swap_double( &(r->x1), &(r->x2) );
      swap_double( &(r->y1), &(r->y2) );
    } 
  
  /* set delimiting angles */
  r->ang_start = atan2( r->y1 - r->cy, r->x1 - r->cx );
  if( r->ang_start < 0 ) r->ang_start += M_2__PI; 
  r->ang_end = atan2( r->y2 - r->cy, r->x2 - r->cx );
  if( r->ang_end < 0 ) r->ang_end += M_2__PI;

  /* if the ends of the last rectangles got interwined, swap the ends 
     of the ring  */
  if( r->ang_start > r->ang_end ) tang = r->ang_end + M_2__PI;
  else tang = r->ang_end;      
  if( (spir > M_PI) && (tang - r->ang_start < M_1_2_PI/2.0) ) 
    swap_ring( r );
  
  /* set ring widths */
  if( conic == CIRCLE ) 
    circ_ring_width( reg, reg_size, r );
  else
    ell_ring_width( reg, reg_size, r );

  /* abnormal case: width of the ring is bigger than little axis -> 
     discard ring */  
  if( r->wmax - r->wmin > r->bx ) return FALSE;
  if( r->ax > msize ) return FALSE;
  
  return TRUE;
}
