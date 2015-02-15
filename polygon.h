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

  polygon.h - This file belongs to ELSDc project (Ellipse and Line Segment 
              Detector with continuous validation).
            - It defines data structures to handle a connected region that 
              can be approximated by a polygon, and prototypes for functions 
              that handle polygons.

------------------------------------------------------------------------------*/


#ifndef POLYGON_H
#define POLYGON_H

#include "misc.h"
#include "rectangle.h"

/*----------------------------------------------------------------------------*/
/** Convex polygon = set of (not necessarily meaningful) rectangles.
 */
typedef struct 
{
  int dim_max;            /* number of maximum allocated number of segments  */
  int dim;                /* number of segments in the polygon               */
  double wmin;            /* dmin is the minimum of all the rectangles' dmin */
  double wmax;            /* dmax is the maximum of all the rectangles' dmax */
  Rectangle* rectlist;    /* list of rectangles composing the polygon        */
} PolyRect;

PolyRect* new_polyrect(void);
void clear_polyrect( PolyRect *poly );
void add_rect_to_polyrect( PolyRect *poly, Rectangle *r );
void write_polyrect( FILE *f, PolyRect *poly );
/* double smooth_score( PolyRect *poly, int *idx, int order, int dim ); */


/*----------------------------------------------------------------------------*/
/** Polygon defined through consecutive ends of its segments.
 */
typedef struct 
{
  int dim;                /* number of segment endpoints in the polygon; it is 
                             twice the number of line segments of the polygon */
  PointD *pts;            /* endpoints of the segments in the polygon         */
} Polygon;

void polyrect2polygon( PolyRect *poly, Polygon *p );

#endif
