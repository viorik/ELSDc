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

  gauss.h - This file belongs to ELSDc project (Ellipse and Line Segment 
            Detector with continuous validation).
          - It defines data structure for Gaussian filter and function 
            prototype for applying it. 

------------------------------------------------------------------------------*/


#ifndef GAUSS_H
#define GAUSS_H


/*----------------------------------------------------------------------------*/
/** Structure for defining the parameters of a Gaussian filter: dimension, mean, 
    standard deviation, and actual values.
 */ 
typedef struct GaussFilter
{
  int dim;
  double sigma;
  double mean;
  double *values;
} *PGaussFilter; /* pointer to Gauss filter structure */

PImageDouble gaussian_sampler( PImageDouble in, double scale, 
                               double sigma_scale );

#endif
