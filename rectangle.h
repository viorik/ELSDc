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

  rectangle.h - This file belongs to ELSDc project (Ellipse and Line Segment 
                Detector with continuous validation).
              - It defines data structures to handle a connected region that 
                can be approximated by a rectangle, and prototypes for functions 
                that handle rectangles (e.g. computing rectangle's parameters).

------------------------------------------------------------------------------*/

#ifndef RECTANGLE_H
#define RECTANGLE_H

/*----------------------------------------------------------------------------*/
/** Rectangle structure. A rectangle is actually a line segment with width.
 */
typedef struct 
{
  double x1, y1, x2, y2;  /* first and second point of the line segment */
  double len;             /* length of the rectangle                    */
  double width;           /* rectangle width                            */
  double wmin, wmax;      /* width towards interior and exterior        */
  double cx, cy;          /* center of the rectangle                    */
  double theta;           /* angle between rectangle orientation and Ox */
  double dx, dy;          /* components of the line segment angle       */
  double prec;            /* tolerance angle                            */
} Rectangle;

void copy_rect( Rectangle *in, Rectangle *out );
void write_rectangle( FILE *f, Rectangle *r );

#endif

