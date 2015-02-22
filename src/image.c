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

  image.c - This file belongs to ELSDc project (Ellipse and Line Segment 
            Detector with continuous validation).
          - It contains functions to free / allocate / allocate and initialize 
            images of type char / int / double . 

------------------------------------------------------------------------------*/


#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include "misc.h"
#include "image.h"


/*----------------------------------------------------------------------------*/
/** Free memory used in PImageChar 'i'.
 */
void free_PImageChar( PImageChar i )
{
  /* check parameters */
  if( (i == NULL) || (i->data == NULL) )
    error("free_PImageChar: invalid input image.");
  free( (void *) i->data );
  free( (void *) i );
}


/*----------------------------------------------------------------------------*/
/** Create a new PImageChar of size 'xsize' times 'ysize'.
 */
PImageChar new_PImageChar( unsigned int xsize, unsigned int ysize )
{
  PImageChar image;

  /* check parameters */
  if( (xsize == 0) || (ysize == 0) ) 
    error("new_PImageChar: invalid image size.");

  /* get memory */
  image = (PImageChar) malloc( sizeof( struct ImageChar ) );
  if( image == NULL ) error("new_PImageChar: not enough memory.");
  image->data = (unsigned char *) calloc( (size_t) (xsize*ysize),
                                          sizeof( unsigned char ) );
  if( image->data == NULL ) error("new_PImageChar: not enough memory.");

  /* set image size */
  image->xsize = xsize;
  image->ysize = ysize;

  return image;
}


/*----------------------------------------------------------------------------*/
/** Create a new PImageChar of size 'xsize' times 'ysize',
    initialized to the value 'fill_value'.
 */
PImageChar new_PImageChar_ini( unsigned int xsize, unsigned int ysize,
                               unsigned char fill_value )
{
  PImageChar image; 
  unsigned int N; 
  unsigned int i;

  /* create image */
  image = new_PImageChar( xsize, ysize ); 

  /* initialize */
  N = xsize * ysize;
  for( i=0; i<N; i++ ) image->data[i] = fill_value;

  return image;
}


/*----------------------------------------------------------------------------*/
/** Free memory used in PImageInt 'i'.
 */
void free_PImageInt( PImageInt i )
{
  /* check parameters */
  if( (i == NULL) || (i->data == NULL) )
    error("free_PImageInt: invalid input image.");
  free( (void *) i->data );
  free( (void *) i );
}


/*----------------------------------------------------------------------------*/
/** Create a new PImageInt of size 'xsize' times 'ysize'.
 */
PImageInt new_PImageInt( unsigned int xsize, unsigned int ysize)
{
  PImageInt image;

  /* check parameters */
  if( (xsize == 0) || (ysize == 0) ) 
    error("new_PImageInt: invalid image size.");

  /* get memory */
  image = (PImageInt) malloc( sizeof( struct ImageInt ) );
  if( image == NULL ) error("new_PImageInt: not enough memory.");
  image->data = (int *) calloc( (size_t) (xsize*ysize), sizeof( int ) );
  if( image->data == NULL ) error("new_PImageInt: not enough memory.");

  /* set image size */
  image->xsize = xsize;
  image->ysize = ysize;

  return image;
}


/*----------------------------------------------------------------------------*/
/** Create a new PImageInt of size 'xsize' times 'ysize',
    initialized to the value 'fill_value'.
 */
PImageInt new_PImageInt_ini( unsigned int xsize, unsigned int ysize,
                             int fill_value )
{
  PImageInt image;
  unsigned int N;
  unsigned int i;

  /* create image */
  image = new_PImageInt( xsize, ysize ); 

  /* initialize */
  N = xsize * ysize;
  for( i=0; i<N; i++ ) image->data[i] = fill_value;

  return image;
}


/*----------------------------------------------------------------------------*/
/** Free memory used in PImageDouble 'i'.
 */
void free_PImageDouble( PImageDouble i )
{
  /* check parameters */
  if( (i == NULL) || (i->data == NULL) )
    error("free_PImageDouble: invalid input image.");
  free( (void *) i->data );
  free( (void *) i );
}


/*----------------------------------------------------------------------------*/
/** Create a new PImageDouble of size 'xsize' times 'ysize'.
 */
PImageDouble new_PImageDouble( unsigned int xsize, unsigned int ysize )
{
  PImageDouble image;

  /* check parameters */
  if( (xsize == 0) || (ysize == 0) ) 
    error("new_PImageDouble: invalid image size.");

  /* get memory */
  image = (PImageDouble) malloc( sizeof( struct ImageDouble ) );
  if( image == NULL ) error("new_PImageDouble: not enough memory.");
  image->data = (double *) calloc( (size_t) (xsize*ysize), sizeof( double ) );
  if( image->data == NULL ) error("new_PImageDouble: not enough memory.");

  /* set image size */
  image->xsize = xsize;
  image->ysize = ysize;

  return image;
}


/*----------------------------------------------------------------------------*/
/** Create a new PImageDouble of size 'xsize' times 'ysize',
    initialized to the value 'fill_value'.
 */
PImageDouble new_PImageDouble_ini( unsigned int xsize, unsigned int ysize,
                                   double fill_value )
{
  PImageDouble image;
  unsigned int N;
  unsigned int i;

  /* create image */
  image = new_PImageDouble( xsize, ysize ); 

  /* initialize */
  N = xsize * ysize;
  for( i=0; i<N; i++ ) image->data[i] = fill_value;

  return image;
}


/*----------------------------------------------------------------------------*/
/** Compute image gradient orientations.
 */
PImageDouble img_gradient_angle( PImageDouble in, double threshold )
{
  unsigned int xsize, ysize, adr, i, j;
  double com1, com2, gx, gy, norm, norm2;
  PImageDouble angles;

  /* check parameters */
  if( in == NULL || in->data == NULL || in->xsize <= 0 || in->ysize <= 0 )
    error("img_gradient_angle: invalid input image");
  if( threshold < 0.0 ) error("img_gradient_angle: threshold must be positive");
  
  xsize = in->xsize;
  ysize = in->ysize;

  /* allocate output image */
  angles = new_PImageDouble(xsize,ysize);
 
  /* 'undefined' on the down and right boundaries */
  for( i=0; i<ysize; i++) angles->data[ i * xsize + xsize - 1 ] = NOTDEF;
  for( j=0; j<xsize; j++) angles->data[ (ysize-1) * xsize + j ] = NOTDEF;

  /*** remaining part ***/
  for( i=0; i<ysize-1; i++ )
    for( j=0; j<xsize-1; j++ )
      {
        adr = i*xsize+j;
        /*
           Norm 2 computation using 2x2 pixel window:
             A B
             C D
           and
             com1 = D-A,  com2 = B-C.
           Then
             gx = B+D - (A+C)   horizontal difference
             gy = C+D - (A+B)   vertical difference
           com1 and com2 are just to avoid 2 additions.
         */
        com1 = in->data[ adr + xsize + 1 ] - in->data[ adr ];
        com2 = in->data[ adr + 1 ] - in->data[ adr + xsize ];
        gx = com1 + com2;
        gy = com1 - com2;
        norm2 = gx*gx+gy*gy;
        norm = sqrt( norm2 / 4.0 );

        if( norm <= threshold ) /* norm too small, gradient not defined */
          angles->data[ adr ] = NOTDEF;
        else    
	  angles->data[ adr ] = atan2( gy, gx );  
      }

  return angles;
}


/*----------------------------------------------------------------------------*/
/** Compute image gradient (magnitude and orientation) and return also the list
    of points ordered by gradient magnitude.
 */
void img_gradient_sort( PImageDouble in, double threshold, CoordList **list_p, 
                        void **mem_p, unsigned int n_bins, double max_grad, 
                        PImageDouble *angles, PImageDouble *gradmag, 
                        PImageDouble *gradx, PImageDouble *grady )
{
  unsigned int xsize, ysize, adr, ind, i, j;
  double com1, com2, gx, gy, norm, norm2;

  /* the rest of the variables are used for pseudo-ordering
     the gradient magnitude values */
  int list_count = 0;
  CoordList *list;
  CoordList **range_l_s; /* array of pointers to start of bin list */
  CoordList **range_l_e; /* array of pointers to end of bin list */
  CoordList *start;
  CoordList *end;

  /* check parameters */
  if( (in == NULL)     || (in->data == NULL) || 
      (in->xsize) <= 0 || (in->ysize <= 0  ) )
    error("img_gradient_sort: invalid input image.");
  if( threshold < 0.0 ) 
    error("img_gradient_sort: 'threshold' must be positive.");
  if( list_p == NULL ) error("img_gradient_sort: NULL pointer 'list_p'.");
  if( mem_p == NULL ) error("img_gradient_sort: NULL pointer 'mem_p'.");
  if( n_bins <= 0 ) error("img_gradient_sort: 'n_bins' must be positive.");
  if( max_grad <= 0.0 ) 
    error("img_gradient_sort: 'max_grad' must be positive.");

  xsize = in->xsize;
  ysize = in->ysize;

  /* allocate output angles image */
  *angles = new_PImageDouble( in->xsize, in->ysize );
  /* allocate image of gradient modulus and gradient vector components
     on Ox and Oy */
  *gradmag = new_PImageDouble( in->xsize, in->ysize );
  *gradx   = new_PImageDouble( in->xsize, in->ysize );
  *grady   = new_PImageDouble( in->xsize, in->ysize );

  /* get memory for "ordered" coordinate list */
  list = (CoordList *) calloc( (size_t) (xsize*ysize), sizeof(CoordList) );
  *mem_p = (void *) list;
  range_l_s = (CoordList **) calloc( n_bins, sizeof(CoordList *) );
  range_l_e = (CoordList **) calloc( n_bins, sizeof(CoordList *) );
  if( (list == NULL) || (range_l_s == NULL) || (range_l_e == NULL) )
    error("ll_angle: not enough memory for sorted list.");

  for( i=0; i<n_bins; i++ ) 
    range_l_s[i] = range_l_e[i] = NULL;

  /* 'undefined' on the down and right boundaries */
  for( i=0; i<ysize; i++ ) (*angles)->data[ i * xsize + xsize - 1 ] = NOTDEF;
  for( j=0; j<xsize; j++ ) (*angles)->data[ (ysize-1) * xsize + j ] = NOTDEF;

  /*** remaining part ***/
  for( i=0; i<ysize-1; i++ )
    for( j=0; j<xsize-1; j++ )
      {
        adr = i * xsize + j;
        /*
           Norm 2 computation using 2x2 pixel window:
             A B
             C D
           and
             com1 = D-A,  com2 = B-C.
           Then
             gx = B+D - (A+C)   horizontal difference
             gy = C+D - (A+B)   vertical difference
           com1 and com2 are just to avoid 2 additions.
         */
        com1 = in->data[ adr + xsize + 1 ] - in->data[ adr ];
        com2 = in->data[ adr + 1 ] - in->data[ adr + xsize ];
        gx = com1 + com2;
        gy = com1 - com2;
        norm2 = gx*gx + gy*gy;
        norm  = sqrt( norm2 / 4.0 );

        (*gradmag)->data[adr] = norm;
	(*gradx)->data  [adr] = gx/2.0;
	(*grady)->data  [adr] = gy/2.0;
        if( norm <= threshold ) /* norm too small, gradient not defined */
          {
            (*angles)->data[adr] = NOTDEF;
          }
        else
	  {
            /* angle computation */
	    (*angles)->data[adr] = atan2( gy, gx );
            /* store the point in the right bin according to its norm */
            ind = (unsigned int) (norm * (double) n_bins / max_grad);
            if( ind >= n_bins ) ind = n_bins - 1;
            if( range_l_e[ind] == NULL)
                range_l_s[ind] = range_l_e[ind] = list + list_count++;
            else
	      {
                range_l_e[ind]->next = list + list_count;
                range_l_e[ind]       = list + list_count++;
              }
            range_l_e[ind]->x    = (int) j;
            range_l_e[ind]->y    = (int) i;
            range_l_e[ind]->next = NULL;
         }
      }

  /* Make the list of points "ordered" by norm value.
     It starts by the larger bin, so the list starts by the
     pixels with higher gradient value. */
  for( i=n_bins-1; (i>0) && (range_l_s[i]==NULL); i-- ) {};
  start = range_l_s[i];
  end   = range_l_e[i];
  if( start != NULL )
    for( i--; i>0 ; i-- )
      if( range_l_s[i] != NULL )
        {
          end->next = range_l_s[i];
          end       = range_l_e[i];
        }
  *list_p = start;

  /* free memory */
  free( (void *) range_l_s );
  free( (void *) range_l_e );
}


/*----------------------------------------------------------------------------*/
/** Mark with the given label the image pixels contained in 'pts' between 
    'start' and 'end' index.
 */
void mark_img_pts( PImageInt im, Point *pts, int start, int end, int label )
{
  int i = 0;
  /* check parameters */
  if( (im == NULL) || (im->data == NULL) || (im->xsize <= 0) || (im->ysize<=0) )
    error("mark_img_pts: invalid input image.");
  if( pts == NULL ) error("mark_img_pts: invalid input point list.");
  if( start >= end ) error("mark_img_pts: invalid index.");
 
  for( i=start; i<end; i++ ) 
    { 
      im->data[ pts[i].y * im->xsize + pts[i].x ] = label; 
    }
}


/*----------------------------------------------------------------------------*/
/** Check if a given point is inside an image.
 */
int in_image( int x, int y, unsigned int xsize, unsigned int ysize)
{
  /* check parameters */
  if( (xsize <= 0) || (ysize<=0) )
    error("in_image: invalid input image dimension.");
  if( (x >= 0) && (x < (int)xsize) && (y >= 0) && (y < (int)ysize) )
    return TRUE;
  return FALSE;
}


/*----------------------------------------------------------------------------*/
/** Check if a given point is a local maximum in an image. Consider 
    8-neighbourhood.
 */
int local_max( PImageDouble im, int x, int y )
{
  int xx, yy;
  int adr = 0;
  /* check parameters */
  if( (im == NULL) || (im->data == NULL) || (im->xsize <= 0) || (im->ysize<=0) )
    error("local_max: invalid input image.");
  adr = y*im->xsize+x;
  for( xx=x-1; xx<=x+1; xx++ )
    for( yy=y-1; yy<=y+1; yy++ )
      if( (xx!=x) && (yy!=y) && (in_image( xx, yy, im->xsize, im->ysize )) )
        if( im->data[yy*im->xsize+xx] > im->data[adr] ) return FALSE;

  return TRUE;
}
