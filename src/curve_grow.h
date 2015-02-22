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


#ifndef CURVE_GROW_H
#define CURVE_GROW_H

void region2rect( PImageDouble gradmag, Point *reg, int start, int end,
                         double reg_angle, Rectangle *rec, double prec );
int curve_grow( PImageDouble gradmag, PImageDouble angles, 
                       PImageInt used, Point *reg, int *reg_size, 
                       double density_th, double prec, PolyRect *poly, 
                       int *label, double *pext1, double *pext2, double *spir );
#endif
