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


  rectangle.c - This file belongs to ELSDc project (Ellipse and Line Segment 
                Detector with continuous validation).
              - It contains functions to work with rectangular connected 
                regions of pixels.

------------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "misc.h"
#include "rectangle.h"

/*----------------------------------------------------------------------------*/
/** Copy one Rectangle structure to another.
 */
void copy_rect( Rectangle *in, Rectangle *out )
{
  if( (in == NULL) || (out == NULL) ) 
    error("copy_rect: invalid 'in' or 'out'.");
  out->x1 = in->x1;
  out->y1 = in->y1;
  out->x2 = in->x2;
  out->y2 = in->y2;
  out->len = in->len;
  out->width = in->width;
  out->wmin = in->wmin;
  out->wmax = in->wmax;
  out->cx = in->cx;
  out->cy = in->cy;
  out->theta = in->theta;
  out->dx = in->dx;
  out->dy = in->dy;
  out->prec = in->prec;
}


/*----------------------------------------------------------------------------*/
/** Write Rectangle coordinates to file; if file is NULL, write to stdout.
 */
void write_rectangle( FILE *f, Rectangle *r )
{
  fprintf(f, "%f %f %f %f %f %f %f %f\n", r->x1, r->y1, r->x2, r->y2, 
                                          r->dx, r->dy, r->wmax, r->wmin ); 
}
