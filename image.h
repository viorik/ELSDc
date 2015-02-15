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

  image.h - This file belongs to ELSDc project (Ellipse and Line Segment 
            Detector with continuous validation).
          - It defines data structures for char / int / double images, function 
            prototypes to free / allocate / allocate and initialize these 
            structures, and prototype for image gradient computation. 

  ----------------------------------------------------------------------------*/

#ifndef IMAGE_H
#define IMAGE_H

#include "misc.h"

/*----------------------------------------------------------------------------*/
/** char image data type.
 */
typedef struct ImageChar
{ 
  unsigned char *data;
  unsigned int xsize, ysize;
} *PImageChar; /* pointer to ImageChar */

void free_PImageChar( PImageChar i );
PImageChar new_PImageChar( unsigned int xsize, unsigned int ysize );
PImageChar new_PImageChar_ini( unsigned int xsize, unsigned int ysize,
                               unsigned char fill_value );


/*----------------------------------------------------------------------------*/
/** int image data type.
 */
typedef struct ImageInt
{
  int *data;
  unsigned int xsize, ysize;
} *PImageInt; /* pointer to ImageInt */

void free_PImageInt( PImageInt i );
PImageInt new_PImageInt( unsigned int xsize, unsigned int ysize );
PImageInt new_PImageInt_ini( unsigned int xsize, unsigned int ysize,
                             int fill_value );


/*----------------------------------------------------------------------------*/
/** double image data type.
 */
typedef struct ImageDouble
{
  double *data;
  unsigned int xsize, ysize;
} *PImageDouble; /* pointer to ImageDouble */


/*----------------------------------------------------------------------------*/
/** Data structure for linked list of point coordinates. 
 */
typedef struct list CoordList;
struct list 
{
  int x, y;
  CoordList *next;
};


void free_PImageDouble( PImageDouble i );
PImageDouble new_PImageDouble( unsigned int xsize, unsigned int ysize );
PImageDouble new_PImageDouble_ini( unsigned int xsize, unsigned int ysize,
                                   double fill_value );

PImageDouble img_gradient_angle( PImageDouble in, double threshold );
void img_gradient_sort( PImageDouble in, double threshold, CoordList **list_p, 
                        void **mem_p, unsigned int n_bins, double max_grad, 
                        PImageDouble *angles, PImageDouble *gradmag, 
                        PImageDouble *gradx, PImageDouble *grady );
void mark_img_pts( PImageInt im, Point *pts, int start, int end, int label );
int in_image( int x, int y, unsigned int xsize, unsigned int ysize);
int local_max( PImageDouble im, int x, int y );
#endif
