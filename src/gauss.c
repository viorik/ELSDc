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

  gauss.c - This file belongs to ELSDc project (Ellipse and Line Segment 
            Detector with continuous validation).
          - It contains functions to free, allocate, initialize, and apply a 
            Gaussian filter on an image.

------------------------------------------------------------------------------*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "image.h"
#include "gauss.h"


/*----------------------------------------------------------------------------*/
/** Free memory used in a PGaussFilter 'in'.
 */
static void free_PGaussFilter( PGaussFilter in )
{
  if( (in == NULL) || (in->values == NULL) )
    error("free_PGaussFilter: invalid filter input.");
  free( (void *) in->values );
  free( (void *) in );
}


/*----------------------------------------------------------------------------*/
/** Allocate space for a new PGaussFilter of size 'dim'.
 */
static PGaussFilter new_PGaussFilter( unsigned int dim )
{
  PGaussFilter fil;

  /* check parameters */
  if( dim <= 0 ) error("new_PGaussFilter: invalid filter size");

  fil = (PGaussFilter) malloc( sizeof( struct GaussFilter ) );
  if( fil == NULL ) error("new_PGaussFilter: not enough memory");
  fil->values = (double *) calloc( dim, sizeof(double) );
  if( fil->values == NULL ) error("new_PGaussFilter: not enough memory");

  fil->dim = dim;
  fil->mean = 0.0;

  return fil;
}


/*----------------------------------------------------------------------------*/
/**
   Compute a Gaussian kernel of length 'kernel->dim',
   standard deviation 'sigma', and centered at value 'mean'.
   For example, if mean=0.5, the Gaussian will be centered
   in the middle point between values 'kernel->values[0]'
   and 'kernel->values[1]'.
 */
static void gaussian_kernel( PGaussFilter kernel, double sigma, double mean )
{
  double sum = 0.0;
  double val;
  int i;

  /* check parameters */
  if( (kernel == NULL) || (kernel->values == NULL) )
    error("gaussian_kernel: invalid struct 'kernel'.");
  if( sigma <= 0.0 ) error("gaussian_kernel: 'sigma' must be positive.");

  /* compute gaussian kernel */
  kernel->mean = mean;
  kernel->sigma = sigma; 

  for( i=0; i<kernel->dim; i++ )
    {
      val = ( (double) i - mean ) / sigma;
      kernel->values[i] = exp( -0.5 * val * val );
      sum += kernel->values[i];
    }

  /* normalization */
  if( sum >= 0.0 ) 
    for( i=0; i<kernel->dim; i++ ) 
      kernel->values[i] /= sum;
}


/*----------------------------------------------------------------------------*/
/**
   Subsample image 'in' with Gaussian filtering, to a scale 'scale'
   (for example, 0.8 will give a result at 80% of the original size),
   using a standard deviation sigma given by:

     sigma = sigma_scale / scale,   if scale <  1.0
     sigma = sigma_scale,           if scale >= 1.0
 */
PImageDouble gaussian_sampler( PImageDouble in, double scale, 
                               double sigma_scale )
{
  PImageDouble aux;
  PImageDouble out;
  PGaussFilter kernel;
  unsigned int N, M, h, n, x, y;
  int i;
  int xc, yc, j, double_x_size, double_y_size;
  double sigma, xx, yy, sum, prec;
  
  /* check parameters */
  if( (in == NULL) || (in->data == NULL) || (in->xsize <= 0) || 
      (in->ysize <= 0) )
    error("gaussian_sampler: invalid image.");
  if( scale <= 0.0 ) error("gaussian_sampler: 'scale' must be positive.");
  if( sigma_scale <= 0.0 )
    error("gaussian_sampler: 'sigma_scale' must be positive.");

  /* get memory for images */
  N = (unsigned int) floor( in->xsize * scale );
  M = (unsigned int) floor( in->ysize * scale );
  aux = new_PImageDouble( N, in->ysize );
  out = new_PImageDouble( N, M );

  /* sigma, kernel size and memory for the kernel */
  sigma = scale < 1.0 ? sigma_scale / scale : sigma_scale;

  /*
     The size of the kernel is selected to guarantee that the
     the first discarded term is at least 10^prec times smaller
     than the central value. For that, h should be larger than x, with
       e^(-x^2/2sigma^2) = 1/10^prec.
     Then,
       x = sigma * sqrt( 2 * prec * ln(10) ).
   */
  prec = 3.0;
  h = (unsigned int) ceil( sigma * sqrt( 2.0 * prec * log(10.0) ) );
  n = 1 + 2*h; /* kernel size */
  kernel = new_PGaussFilter(n);

  /* auxiliary double image size variables */
  double_x_size = (int) (2 * in->xsize);
  double_y_size = (int) (2 * in->ysize);

  /* First subsampling: x axis */
  for( x=0; x<aux->xsize; x++ ) 
    {
      /*
         x   is the coordinate in the new image.
         xx  is the corresponding x-value in the original size image.
         xc  is the integer value, the pixel coordinate of xx.
       */
      xx = (double) x / scale;
      /* coordinate (0.0,0.0) is in the center of pixel (0,0),
         so the pixel with xc=0 get the values of xx from -0.5 to 0.5 */
      xc = (int) floor( xx + 0.5 );
      gaussian_kernel( kernel, sigma, (double) h + xx - (double) xc );
      /* the kernel must be computed for each x because the fine
         offset xx-xc is different in each case */

      for( y=0; y<aux->ysize; y++ )
        {
          sum = 0.0;
          for( i=0; i<kernel->dim; i++ )
            {
              j = xc - h + i;

              /* symmetry boundary condition */
              while( j < 0 ) j += double_x_size;
              while( j >= double_x_size ) j -= double_x_size;
              if( j >= (int) in->xsize ) j = double_x_size-1-j;

              sum += in->data[ j + y * in->xsize ] * kernel->values[i];
            }
          aux->data[ x + y * aux->xsize ] = sum;
        }
    }

  /* Second subsampling: y axis */
  for( y=0; y<out->ysize; y++ )
    {
      /*
         y   is the coordinate in the new image.
         yy  is the corresponding x-value in the original size image.
         yc  is the integer value, the pixel coordinate of xx.
       */
      yy = (double) y / scale;
      /* coordinate (0.0,0.0) is in the center of pixel (0,0),
         so the pixel with yc=0 get the values of yy from -0.5 to 0.5 */
      yc = (int) floor( yy + 0.5 );
      gaussian_kernel( kernel, sigma, (double) h + yy - (double) yc );
      /* the kernel must be computed for each y because the fine
         offset yy-yc is different in each case */

      for( x=0; x<out->xsize; x++ )
        {
          sum = 0.0;
          for( i=0; i<kernel->dim; i++ )
            {
              j = yc - h + i;

              /* symmetry boundary condition */
              while( j < 0 ) j += double_y_size;
              while( j >= double_y_size ) j -= double_y_size;
              if( j >= (int) in->ysize ) j = double_y_size-1-j;

              sum += aux->data[ x + j * aux->xsize ] * kernel->values[i];
            }
          out->data[ x + y * out->xsize ] = sum;
        }
    }

  /* free memory */
  free_PGaussFilter(kernel);
  free_PImageDouble(aux);
  return out;
}
