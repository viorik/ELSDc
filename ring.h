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

  ring.h - This file belongs to ELSDc project (Ellipse and Line Segment 
           Detector with continuous validation).
         - It defines data structures to handle a connected region that 
           can be approximated by a circular or elliptical ring, and prototypes 
           for functions that handle rings (e.g. computing ring's parameters, 
           check validity, compute distance between a point and a ring etc).

------------------------------------------------------------------------------*/


#ifndef RING_H
#define RING_H

#define CIRCLE 1
#define ELLIPSE 0

/*----------------------------------------------------------------------------*/
/** Elliptical/circular ring structure. A ring is defined by end points, width,
    center, orientation, axes, delimiting angles. 
 */
typedef struct 
{
  double x1, y1, x2, y2;     /* End points of the circle/ellipse arc */
  double width;              /* ring width */
  double cx, cy;             /* center of the circle/ellipse */
  double theta;              /* ellipse orientation; 0 for circle */
  double ax, bx;             /* ellipse axes; ax=bx for circle */
  double ang_start, ang_end; /* delimiting angles */
  double wmin, wmax;         /* width towards interior and exterior */ 
  int full;                  /* if full == 1, circle/ellipse complete, 
                                used for display */
} Ring;

void copy_ring( Ring *in, Ring *out );
void swap_ring( Ring *r );
int check_circ_ring( Ring *cring );
int check_ell_ring( Ring *ering );
int is_in_circ_ring( Ring * r, int x, int y );
int is_in_ell_ring( Ring * r, int x, int y );
void rosin_point( Ring *ering, int x, int y, double *xi, double *yi );
double d_rosin( Ring *ering, double x, double y );
void ellipse_foci( double *param, double *foci );
double ellipse_normal_angle( double x, double y, double *foci );
int get_ring( Point *reg, int reg_size, double ang, double *param, 
              double *pext1, double *pext2, int conic, double spir, Ring *r, 
              int *grad_dir, double *foci, int msize );
#endif
