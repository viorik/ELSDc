/*------------------------------------------------------------------------------
  
  Copyright (c) 2007-2011 rafael grompone von gioi (grompone@gmail.com)
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


  iterator.c - This file belongs to ELSDc project (Ellipse and Line Segment 
               Detector with continuous validation).
             - It contains functions to scan the points in a rectangular area.

------------------------------------------------------------------------------*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include "misc.h"
#include "rectangle.h"
#include "iterator.h"


/*----------------------------------------------------------------------------*/
/** Get abscise of the intersection point between lower side of the rectangle 
    and the vertical line.
 */
static double inter_low( double x, double x1, double y1, double x2, double y2 )
{
  if( (x1 > x2) || (x < x1) || (x > x2) )
    {
      fprintf( stderr,"inter_low: x %g x1 %g x2 %g.\n", x, x1, x2 );
      error("inter_low: something went wrong!");
    }
  if( (double_equal( x1, x2 )) && (y1<y2) ) return y1;
  if( (double_equal( x1, x2 )) && (y1>y2) ) return y2;
  return y1 + (x-x1) * (y2-y1) / (x2-x1);
}


/*----------------------------------------------------------------------------*/
/** Get abscise of the intersection point between higher side of the rectangle 
    and the vertical line.
 */
static double inter_hi( double x, double x1, double y1, double x2, double y2 )
{
  if( (x1 > x2) || (x < x1) || (x > x2) )
    {
      fprintf(stderr,"inter_hi: x %g x1 %g x2 %g.\n",x,x1,x2);
      error("inter_hi: something went wrong!");
    }
  if( (double_equal( x1, x2 )) && (y1<y2) ) return y2;
  if( (double_equal( x1, x2 )) && (y1>y2) ) return y1;
  return y1 + (x-x1) * (y2-y1) / (x2-x1);
}


/*----------------------------------------------------------------------------*/
/** Free memory used by a rectangle iterator.
 */
void free_RectIter( RectIter *iter )
{
  /* check parameters */
  if( iter == NULL ) error("free_RectIter: NULL iterator.");
  free( (void *) iter );
}


/*----------------------------------------------------------------------------*/
/** Check if the iterator finished the full iteration.
 */
int end_RectIter( RectIter *i )
{
  /* check parameters */
  if( i == NULL ) error("end_RectIter: NULL iterator.");
  return (double)(i->x) > i->vx[2];
}


/*----------------------------------------------------------------------------*/
/** Increment a rectangle iterator.
 */
void inc_RectIter( RectIter *i )
{
  /* check parameters */
  if( i == NULL ) error("inc_RectIter: NULL iterator.");

  if( (double) (i->x) <= i->vx[2] ) i->y++;
  while( ((double) (i->y) > i->ye) && ((double) (i->x) <= i->vx[2]) )
    {
      /* new x */
      i->x++;

      if( (double) (i->x) > i->vx[2] ) return; /* end of iteration */

      /* update lower y limit for the line */
      if( (double) i->x < i->vx[3] )
        i->ys = inter_low((double)i->x,i->vx[0],i->vy[0],i->vx[3],i->vy[3]);
      else 
        i->ys = inter_low((double)i->x,i->vx[3],i->vy[3],i->vx[2],i->vy[2]);

      /* update upper y limit for the line */
      if( (double)i->x < i->vx[1] )
        i->ye = inter_hi((double)i->x,i->vx[0],i->vy[0],i->vx[1],i->vy[1]);
      else 
        i->ye = inter_hi((double)i->x,i->vx[1],i->vy[1],i->vx[2],i->vy[2]);

      /* new y */
      i->y = (int) ceil(i->ys);
    }
}


/*----------------------------------------------------------------------------*/
/** Create and initialize a rectangle iterator.
 */
RectIter *ini_RectIter( Rectangle *r )
{
  double vx[4], vy[4];
  int n, offset;
  RectIter *i;

  /* check parameters */
  if( r == NULL ) error("ini_RectIter: invalid input rectangle.");

  i = (RectIter *) malloc( sizeof(RectIter) );
  if( i == NULL ) error("ini_RectIter: not enough memory.");

  vx[0] = r->x1 - r->dy * r->wmax;
  vy[0] = r->y1 + r->dx * r->wmax;
  vx[1] = r->x2 - r->dy * r->wmax;
  vy[1] = r->y2 + r->dx * r->wmax;
  vx[2] = r->x2 - r->dy * r->wmin;
  vy[2] = r->y2 + r->dx * r->wmin;
  vx[3] = r->x1 - r->dy * r->wmin;
  vy[3] = r->y1 + r->dx * r->wmin;

  if     ( (r->x1 <  r->x2) && (r->y1 <= r->y2) ) offset = 0;
  else if( (r->x1 >= r->x2) && (r->y1 <  r->y2) ) offset = 1;
  else if( (r->x1 >  r->x2) && (r->y1 >= r->y2) ) offset = 2;
  else offset = 3;

  for( n=0; n<4; n++ )
    {
      i->vx[n] = vx[ (offset+n) % 4 ];
      i->vy[n] = vy[ (offset+n) % 4 ];
    }

  /* starting point */
  i->x = (int) ceil(i->vx[0]) - 1;
  i->y = (int) ceil(i->vy[0]);
  i->ys = i->ye = -DBL_MAX;

  /* advance to the first point */
  inc_RectIter(i);

  return i;
}
