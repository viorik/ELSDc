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

  iterator.h - This file belongs to ELSDc project (Ellipse and Line Segment 
               Detector with continuous validation).
             - It defines a data structure to scan points in a rectangular area,
               and prototypes for functions to use it. 

------------------------------------------------------------------------------*/


#ifndef ITERATOR_H
#define ITERATOR_H

/*----------------------------------------------------------------------------*/
/** Rectangle iterator structure: useful to scan the pixels inside a rectangle.
 */
typedef struct 
{
  double vx[4];
  double vy[4];
  double ys, ye;
  int    x, y;
} RectIter;

void free_RectIter( RectIter *iter );
int end_RectIter( RectIter *i );
void inc_RectIter( RectIter *i );
RectIter *ini_RectIter( Rectangle *r );

#endif
