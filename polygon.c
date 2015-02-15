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


  polygon.c - This file belongs to ELSDc project (Ellipse and Line Segment 
              Detector with continuous validation).
            - It contains functions to work with rectangular connected 
              regions of pixels.

------------------------------------------------------------------------------*/


#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include "misc.h"
#include "rectangle.h"
#include "polygon.h"

/*----------------------------------------------------------------------------*/
/** Add new rectangle to polygon.
 */
PolyRect* new_polyrect( void )
{
  PolyRect* poly;
  poly = (PolyRect*) malloc ( sizeof(PolyRect) );
  poly->dim = 0;
  poly->dim_max = 0;
  poly->wmin = DBL_MAX;
  poly->wmax = -DBL_MAX;
  poly->rectlist = NULL;
  return poly;
}


/*----------------------------------------------------------------------------*/
/** Reset a given polygon.
 */
void clear_polyrect( PolyRect *poly )
{
  /* check parameters */
  if( poly == NULL ) error("clear_poly: polygon should be already malloc'ed.");
 
  poly->dim = 0;
  poly->dim_max = 0;
  poly->wmin = DBL_MAX;
  poly->wmax = -DBL_MAX;
  free(poly->rectlist);
  poly->rectlist = NULL;  
}


/*----------------------------------------------------------------------------*/
/** Add new rectangle to polygon.
 */
void add_rect_to_polyrect( PolyRect *poly, Rectangle *r )
{
  /* if the pre-allocated size of the polygon is overpassed, double the size 
     of the allocated region */
  if( poly->dim > poly->dim_max )
    { 
      poly->rectlist = (Rectangle *) realloc ( poly->rectlist, (2*poly->dim) *
                                               sizeof(Rectangle) );
      poly->dim_max *= 2;
    }
  /* add rectangle to polygon */
  copy_rect( r, &(poly->rectlist[poly->dim-1]) );

  /* update polygon's widths wmin and wmax if new rectangle's wmin and wmax 
     exceed them */
  if( poly->wmin > r->wmin ) poly->wmin = r->wmin;
  if( poly->wmax < r->wmax ) poly->wmax = r->wmax;
}


/*----------------------------------------------------------------------------*/
/** Write polygon to file; if file is NULL, write to stdout.
 */
void write_polyrect( FILE *f, PolyRect *poly )
{
  int i;
  for( i=0; i<poly->dim; i++ )
    write_rectangle( f, &(poly->rectlist[i]) );
}


/*----------------------------------------------------------------------------*/
/** Detect continuity breaks (e.g. corners) in tangent space: represent each 
    segment of the polygon as a horizontal segment in a (length, angle) space, 
    and verify the centres of 3 consecutive rectangles in the polygon for 
    approximate colinearity. 
    'idx': indexes of the rectangles to check;
    'order': 1 if consecutive rectangles, 0 otherwise; 
    'dim' is the size of 'idx';
    Output: continuity score in [0,1]: 0 not smooth, 1 perfect smooth. 
 */
/*double smooth_score( PolyRect *poly, int *idx, int order, int dim )
{
  if( poly == NULL ) 
    error("smooth_score: invalid input."); 
  if( idx == NULL ) 
    error("smooth_score: invalid index list.");
  if( dim<=0 ) error("smooth_score: invalid dim.");
  
  if( dim < 3 ) return 1.0;

  PointD p[3];
  Rectangle r0, r1, r2;
  r0 = poly->rectlist[idx[0]];
  r1 = poly->rectlist[idx[1]];
  r2 = poly->rectlist[idx[2]];
  p[0].x = 0; 
  p[0].y = r0.theta;
  p[1].x = r0.len/2.0 + r1.len/2.0; 
  p[1].y = r1.theta;
  if( order ) p[2].x = r0.len/2.0 + r1.len + r2.len/2.0 ;
  else p[2].x = -r0.len/2.0 - r2.len/2.0 ;
  p[2].y = r2.theta;
}*/


/*----------------------------------------------------------------------------*/
/** Given a polygon as a list of rectangles, extract only the ending points.
    A polygon with 'dim' segments will have '2*dim' ending points.  
 */
void polyrect2polygon( PolyRect *poly, Polygon *p )
{
  int i;
  /* check parameters */
  if( poly == NULL || p == NULL ) 
    error("polyrect2polygon: invalid input.");

  /* allocate memory for polygon */
  p->dim = 2 * poly->dim;
  p->pts = (PointD *) malloc ( p->dim * sizeof(PointD) );
  if( p->pts == NULL ) 
    error("polyrect2polygon: not enough memory.");

  /* copy ending points */
  for( i=0; i<poly->dim; i++ )
    {
      p->pts[2*i  ].x = poly->rectlist[i].x1 + 0.5;
      p->pts[2*i  ].y = poly->rectlist[i].y1 + 0.5;
      p->pts[2*i+1].x = poly->rectlist[i].x2 + 0.5;
      p->pts[2*i+1].y = poly->rectlist[i].y2 + 0.5;
    }
}
