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


  elsdc.c - This file belongs to ELSDc project (Ellipse and Line Segment 
            Detector with continuous validation).
          - It contains the main detection procedure: curve growing and curve
            validation. 

------------------------------------------------------------------------------*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include "misc.h"
#include "image.h"
#include "gauss.h"
#include "rectangle.h"
#include "iterator.h"
#include "polygon.h"
#include "ring.h"
#include "ellipse_fit.h"
#include "curve_grow.h"
#include "elsdc.h"


/*----------------------------------------------------------------------------*/
/** Continuous formulation of the number of false alarms 
    NFAC = NT * k^n / n!
    log(n!) is bounded by Stirling's approximation:
       n! >= sqrt(2pi) * n^(n+0.5) * exp(-n)
     then, log10(NFA) <= log10(NT) + n*log10(k) - log10(latter expansion) 
 */ 
#define LOG10_1_2_M_2__PI  0.39908993418
#define LOG10_EXP1         0.43429448190
static double NFAc( int n, double k, double logNT )
{
  /* check parameters */
  if( n < 0) error("NFAc: 'n' must be strictly positive.");
  if( n==0 ) return -logNT;
  double logNFAC = logNT;
  
  k = max(k,0.01);
  logNFAC = logNT + n * log10(k) - LOG10_1_2_M_2__PI - 
            (n+0.5) * log10(n) + n * LOG10_EXP1 ;

  return -logNFAC;
}


/*----------------------------------------------------------------------------*/
/*----------------- FUNCTIONS FOR RECTANGLE VALIDATION -----------------------*/
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/** Compute additive normalised per-pixel error inside a rectangle, and stack
    visited points so they can be marked as 'used' if the rectangle is the most 
    meaningful primitive.
 */
static void rect_count_error( PImageDouble angles, PImageInt used, 
                              Rectangle *rec, int *pts, double *err, 
                              Point *buff, int *size_buff )
{
  RectIter *i;
  int adr;
  double theta;
 
  /* check parameters */
  if( (angles == NULL) || (angles->data == NULL) )
    error("rect_count_error: invalid input 'angles' image.");
  if( (used == NULL) || (used->data == NULL) )
    error("rect_count_error: invalid input 'used' image.");
  if( rec == NULL ) error("rect_count_error: invalid input rectangle.");
  if( buff == NULL ) error("rect_count_error: invalid buffer.");
  if( size_buff == NULL ) 
    error("rect_count_error: size_buff must be non null.");

  theta = rec->theta - M_1_2_PI;
  *pts = 0; *err = 0.0;

  /* iterate through rectangle and accumulate per-pixel normalised error */
  for( i=ini_RectIter(rec); !end_RectIter(i); inc_RectIter(i) )
    if( in_image( i->x, i->y, used->xsize, used->ysize ) ) 
      {
        ++(*pts);
        adr = (i->y) * used->xsize + (i->x);
        if( (used->data[adr] == USED) || (angles->data[adr] == NOTDEF) )
          (*err) += 1.0;
        else
          (*err) += norm_angle_diff( angles->data[adr], theta );
        /* store points in buffer */
        buff[*size_buff].x = i->x;
        buff[*size_buff].y = i->y;
        ++(*size_buff);
      }
  free_RectIter(i);
}


/*----------------------------------------------------------------------------*/
/** Given a list of points, estimate segment. One segment is a particular case 
    of a polygon.
 */
static void get_seg( PImageDouble angles, PImageDouble gradmag, 
                     PImageInt used, PolyRect *poly, Point *reg, 
                     int reg_size, double prec )
{
  double sumdx = 0.0, sumdy = 0.0;
  double reg_angle;
  int i;
  int adr;
  Rectangle rec;
  /* check parameters */
  if( (angles == NULL) || (angles->data == NULL) )
    error("get_seg: invalid input 'angles' image.");
  if( (gradmag == NULL) || (gradmag->data == NULL) )
    error("get_seg: invalid input 'angles' image.");
  if( (used == NULL) || (used->data == NULL) )
    error("get_seg: invalid input 'used' image.");
  if( poly == NULL ) error("get_seg: invalid input polygon.");
  if( reg == NULL ) error("get_seg: invalid input region.");
  if( reg_size <= 0 ) error("get_seg: 'reg_size' must be strictly positive.");
  
  /* compute the mean angle of the whole region */
  for( i=0; i<reg_size; i++ )
    {
      adr = reg[i].x + reg[i].y * angles->xsize;
      sumdx += cos( angles->data[adr] );
      sumdy += sin( angles->data[adr] );
    }
  reg_angle = atan2( sumdy, sumdx );

  /* estimate rectangle */
  region2rect( gradmag, reg, 0, reg_size, reg_angle, &rec, prec ); 
  
  /* add rectangle to polygon */
  poly->dim++;
  add_rect_to_polyrect( poly, &rec );  
}


/*----------------------------------------------------------------------------*/
/*------------------- FUNCTIONS FOR POLYGON VALIDATION -----------------------*/
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/** Compute additive normalised per-pixel error inside a polygon, and stack
    visited points so they can be marked as 'used' if the polygon is the most 
    meaningful primitive.
 */
static void poly_count_error( PImageDouble angles, PImageInt used, 
                              PolyRect *poly, double polywmin, double polywmax,
                              int *pts, double *err, Point *buff, 
                              int *size_buff ) 
{
  int pts_rect;
  double err_rect;
  int i;
  /* check parameters */
  if( (angles == NULL) || (angles->data == NULL) )
    error("poly_count_error: invalid input 'angles' image.");
  if( (used == NULL) || (used->data == NULL) )
    error("poly_count_error: invalid input 'used' image.");
  if( poly == NULL ) error("poly_count_error: invalid input polygon.");
  if( buff == NULL ) error("poly_count_error: invalid buffer.");
  if( size_buff == NULL ) 
    error("poly_count_error: size_buff must be non null.");

  *pts = 0; *err = 0.0;
  *size_buff = 0;

  /* Scan every rectangle in the polygon and count the additive error. 
     We must consider the same widths wmin and wmax for all rectangles. */
  for( i=0; i<poly->dim; i++)
    { 
      poly->rectlist[i].wmin = polywmin;
      poly->rectlist[i].wmax = polywmax;
      rect_count_error( angles, used, &(poly->rectlist[i]), &pts_rect, 
                        &err_rect, buff, size_buff );
      (*pts) += pts_rect;
      (*err) += err_rect;
    }  
}


/*----------------------------------------------------------------------------*/
/** Iteration to improve NFA: count additive error and compute NFA. 
 */
static double poly_nfa( PImageDouble angles, PImageInt used, PolyRect *poly, 
                        double polywmin, double polywmax, Point *buff, 
                        int *size_buff, double logNT_poly )
{
  int pts;
  double err;
  
  /* check parameters */
  if( (angles == NULL) || (angles->data == NULL) )
    error("poly_nfa: invalid input 'angles' image.");
  if( (used == NULL) || (used->data == NULL) )
    error("poly_nfa: invalid input 'used' image.");
  if( poly == NULL ) error("poly_nfa: invalid input polygon.");
  if( buff == NULL ) error("poly_nfa: invalid buffer.");
  if( size_buff == NULL ) 
    error("poly_nfa: size_buff must be non null.");

  /* count additive error */
  poly_count_error( angles, used, poly, polywmin, polywmax, &pts, &err, buff, 
                    size_buff );
  /* return computed NFA */
  return NFAc( pts, err, logNT_poly );
}


/*----------------------------------------------------------------------------*/
/** Try some polygon variations to improve NFAc value.
 */
static double poly_improve( PImageDouble angles, PImageInt used, PolyRect *poly,
                            Point **local_buff, int *size_local_buff, 
                            Point **tmp_buff, double logNT_poly )
{
  double delta = 0.5;
  double log_nfa,log_nfa_new;

  int size_tmp_buff;
  double wmin,wmax;
  
  /* check parameters */
  if( (angles == NULL) || (angles->data == NULL) )
    error("poly_improve: invalid input 'angles' image.");
  if( (used == NULL) || (used->data == NULL) )
    error("poly_improve: invalid input 'used' image.");
  if( poly == NULL ) error("poly_improve: invalid input polygon.");
  if( *local_buff == NULL ) error("poly_improve: invalid local buffer.");
  if( size_local_buff == NULL ) 
    error("poly_improve: 'size_local_buff' must be non null.");
  if( *tmp_buff == NULL ) error("poly_improve: invalid tmp buffer.");
 
  *size_local_buff = 0;
  wmax = poly->wmax;
  wmin = poly->wmin;
  /* check if polygon valid: wmin and wmax must have opposite sign */
  if ( wmin*wmax >0 ) return -logNT_poly;

  log_nfa = poly_nfa( angles, used, poly, wmin, wmax, *local_buff, 
                      size_local_buff, logNT_poly );
     
  /* try to improve the polygon, by reducing width; first reduce wmax */ 
  /*for( n=0; n<10; n++ )
    {*/      
      while( wmax - delta > 0.0)
        {
          wmax -= delta;
          size_tmp_buff = 0;
  
          log_nfa_new = poly_nfa( angles, used, poly, wmin, wmax, *tmp_buff, 
                                  &size_tmp_buff, logNT_poly );
          if( log_nfa_new > log_nfa )
            {
              log_nfa = log_nfa_new;
              /* switch addresses of the pointers to the local and tmp buffer */
              /*swap_ptr_pts( local_buff, tmp_buff );
              (*size_local_buff) = size_tmp_buff;*/
              poly->wmax -= delta;
            }  
        }
    /* } */

  /* now try to reduce dmin */
  /* for( n=0; n<10; n++ )
    { */      
      while( wmin + delta < 0.0)
        {
          wmin += delta;
          size_tmp_buff = 0;
        
          log_nfa_new = poly_nfa( angles, used, poly, wmin, poly->wmax, 
                                  *tmp_buff, &size_tmp_buff, logNT_poly );
          if( log_nfa_new > log_nfa )
            {
              log_nfa = log_nfa_new;
              /* switch addresses of the pointers to the local and tmp buffer */
              /* swap_ptr_pts( local_buff, tmp_buff );
              (*size_local_buff) = size_tmp_buff; */
              poly->wmin += delta;
            }    
        }
    /* } */
  return log_nfa;
}


/*----------------------------------------------------------------------------*/
/*------------------- FUNCTIONS FOR CIRCLE RING VALIDATION -------------------*/
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/** Compute additive normalised per-pixel error inside a circular ring, and 
    stack visited points so they can be marked as 'used' if the circle ring is
    the most meaningful primitive.
 */
static void circ_count_error( PImageDouble angles, PImageInt used, Ring *cring,
                              int x, int y, int grad_dir, Point *buff, 
                              int *size_buff, int label, int *pts, double *err )
{
  int xx, yy, i, adr;
  double theta;
 
  /* check parameters */
  if( (angles == NULL) || (angles->data == NULL) )
    error("circ_count_error: invalid input 'angles' image.");
  if( (used == NULL) || (used->data == NULL) )
    error("circ_count_error: invalid input 'used' image.");
  if( cring == NULL ) error("circ_count_error: invalid input ring.");
  if( buff == NULL ) error("circ_count_error: invalid local buffer.");
  if( size_buff == NULL ) 
    error("circ_count_error: 'size_local_buff' must be non null.");
  if( label <= 1 ) error("circ_count_error: forbidden label value.");
  
  *pts = 0; *err = 0.0;
  
  /* To compute the error, visit each pixel in the ring in a region growing 
     manner */

  /* First make sure the seed point (x,y) is inside the ring; if not, visit 
     neighbours in a 3x3 neighbourhood to find a point inside the ring */
  *size_buff = 0;
  if( is_in_circ_ring( cring, x, y ) )
    {
      buff[0].x = x;
      buff[0].y = y;
      used->data[ buff[0].y * used->xsize + buff[0].x ] = label;
      *size_buff = 1;
    }
  else
    {
      for( xx=x-2; xx<=x+2; xx++ )
        {
          for( yy=y-2; yy<=y+2; yy++ )
            if( in_image( xx, yy, used->xsize, used->ysize ) && is_in_circ_ring( cring, xx, yy ) )
              {
                buff[0].x = xx;
                buff[0].y = yy;
                used->data[ buff[0].y * used->xsize + buff[0].x ] = label;
                *size_buff = 1; 
                break;
              }
          if( *size_buff == 1 ) break;
        } 
    }
  if( *size_buff == 0 ) return; 

  for( i=0; i < *size_buff; i++ )
    for( xx=buff[i].x-1; xx<=buff[i].x+1; xx++ )
      for( yy=buff[i].y-1; yy<=buff[i].y+1; yy++ )
        if( in_image( xx, yy, used->xsize, used->ysize ) )
          {
            adr = xx + yy * used->xsize;
            /* check if the point has already been used or visited */
            if( (used->data[adr] != USED) && (used->data[adr] != label) ) 
              {
                /* check if point inside the circle ring */
                if( is_in_circ_ring( cring, xx, yy ) )
                  {
                    /* check if gradient converges or diverges when computing 
                       the ideal angle of a point on the circle */
                    theta = atan2( (double)yy - cring->cy, 
                                   (double)xx - cring->cx );
                    if( grad_dir == 0 )
                      {
		        if( theta > 0 ) theta = -(M_PI-theta);
		        else theta = M_PI + theta;			
	              }
                    (*err) += norm_angle_diff( angles->data[adr], theta );
                    used->data[adr] = label;
                    buff[*size_buff].x = xx;
                    buff[*size_buff].y = yy;
                    ++(*size_buff); 
                  }
              }
          } 
  *pts = *size_buff;         
}


/*----------------------------------------------------------------------------*/
/** Iteration to improve NFA: count additive error and compute NFA. 
 */
static double circ_ring_nfa( PImageDouble angles, PImageInt used, Ring *cring,
                             int x, int y, int grad_dir, Point *buff, 
                             int *size_buff, int label, double logNT_ell )
{
  int pts;
  double err;
 
  /* check parameters */
  if( (angles == NULL) || (angles->data == NULL) )
    error("circ_ring_nfa: invalid input 'angles' image.");
  if( (used == NULL) || (used->data == NULL) )
    error("circ_ring_nfa: invalid input 'used' image.");
  if( cring == NULL ) error("circ_ring_nfa: invalid input ring.");
  if( buff == NULL ) error("circ_ring_nfa: invalid local buffer.");
  if( size_buff == NULL ) 
    error("circ_ring_nfa: 'size_local_buff' must be non null.");
  if( label <= 1 ) error("circ_ring_nfa: forbidden label value.");

  /* count additive error */
  circ_count_error( angles, used, cring, x, y, grad_dir, buff, size_buff, 
                    label, &pts, &err );

  /* return computed NFA */
  return NFAc( pts, err, logNT_ell );
}


/*----------------------------------------------------------------------------*/
/** Try some circle ring variations to improve NFAc value.
 */
static double circ_ring_improve( PImageDouble angles, PImageInt used, 
                                 Ring *cring, int x, int y, int grad_dir, 
                                 Point **local_buff, int *size_local_buff, 
                                 Point **tmp_buff, int *label, double logNT_ell)
{    
  double log_nfa, log_nfa_new;
  double delta = 0.5;

  int size_tmp_buff;
  
  /* check parameters */
  if( (angles == NULL) || (angles->data == NULL) )
    error("circ_ring_improve: invalid input 'angles' image.");
  if( (used == NULL) || (used->data == NULL) )
    error("circ_ring_improve: invalid input 'used' image.");
  if( cring == NULL ) error("circ_ring_improve: invalid input ring.");
  if( *local_buff == NULL ) error("circ_ring_improve: invalid local buffer.");
  if( *tmp_buff == NULL ) error("circ_ring_improve: invalid tmp buffer.");
  if( size_local_buff == NULL ) 
    error("circ_ring_improve: 'size_local_buff' must be non null.");
  if( *label <= 1 ) error("circ_ring_improve: forbidden label value.");

  if ( ! check_circ_ring(cring) ) return -logNT_ell;
  (*label)++;
  *size_local_buff = 0;
  log_nfa = circ_ring_nfa( angles, used, cring, x, y, grad_dir, *local_buff, 
                           size_local_buff, *label, logNT_ell );
     
  /* try to improve the ring, by reducing width; first reduce wmax */
  Ring r; 
  copy_ring( cring, &r );
  /*for( n=0; n<10; n++ )
    { */        
      while( r.wmax-delta > 0.0)
        {
          r.wmax -= delta;
          size_tmp_buff = 0;
          (*label)++;
         
          log_nfa_new = circ_ring_nfa( angles, used, &r, x, y, grad_dir,
                                       *tmp_buff, &size_tmp_buff, *label,
                                       logNT_ell );
          if( log_nfa_new > log_nfa )
            {
              cring->wmax -= delta;
              log_nfa = log_nfa_new;
              /* switch addresses of the pointers to the local and tmp buffer */
              /*swap_ptr_pts( local_buff, tmp_buff );
              (*size_local_buff) = size_tmp_buff; */
            }            
        }
    /* } */

  /* now try to reduce wmin */
  copy_ring( cring, &r );
  /*for( n=0; n<10; n++ )
    { */         
      while( r.wmin + delta < 0.0)
        {
          r.wmin += delta;
          size_tmp_buff = 0;
          (*label)++;
         
          log_nfa_new = circ_ring_nfa( angles, used, &r, x, y, grad_dir,
                                       *tmp_buff, &size_tmp_buff, *label,
                                       logNT_ell );
          if( log_nfa_new > log_nfa )
            {
              cring->wmin += delta;
              log_nfa = log_nfa_new;
              /* switch addresses of the pointers to the local and tmp buffer */
              /*swap_ptr_pts( local_buff, tmp_buff );
              (*size_local_buff) = size_tmp_buff;*/
            }            
        }  
    /*} */  
  if ( ! check_circ_ring(cring) ) return -logNT_ell;
  return log_nfa;
}


/*----------------------------------------------------------------------------*/
/*-------------------- FUNCTIONS FOR ELLIPSE RING VALIDATION -----------------*/
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/** Compute additive normalised per-pixel error inside an elliptical ring, and 
    stack visited points so they can be marked as 'used' if the ellipse ring is
    the most meaningful primitive.
 */
static void ell_count_error( PImageDouble angles, PImageInt used, Ring *ering,
                             double *foci, int x, int y, int grad_dir, 
                             Point *buff, int *size_buff, int label, int *pts,
                             double *err )
{
  int xx, yy, i, adr;
  double theta;
  
  /* check parameters */
  if( (angles == NULL) || (angles->data == NULL) )
    error("ell_count_error: invalid input 'angles' image.");
  if( (used == NULL) || (used->data == NULL) )
    error("ell_count_error: invalid input 'used' image.");
  if( ering == NULL ) error("ell_count_error: invalid input ring.");
  if( foci == NULL ) error("ell_count_error: invalid input foci.");
  if( buff == NULL ) error("ell_count_error: invalid local buffer.");
  if( size_buff == NULL ) 
    error("ell_count_error: 'size_local_buff' must be non null.");
  if( label <= 1 ) error("ell_count_error: forbidden label value.");
  
  *pts = 0; *err = 0.0;

  /* To compute the error, visit each pixel in the ring in a region growing 
     manner */
  /* First make sure the seed point (x,y) is inside the ring; if not, visit 
     neighbours in a 3x3 neighbourhood to find a point inside the ring */
  *size_buff = 0;
  if( is_in_ell_ring( ering, x, y ) )
    {
      buff[0].x = x;
      buff[0].y = y;
      used->data[ buff[0].y * used->xsize + buff[0].x ] = label;
      *size_buff = 1;
    }
  else
    {
      for( xx=x-2; xx<=x+2; xx++ )
        {
          for( yy=y-2; yy<=y+2; yy++ )
            if( in_image( xx, yy, used->xsize, used->ysize ) && is_in_ell_ring( ering, xx, yy ) )
              {
                buff[0].x = xx;
                buff[0].y = yy;
                used->data[ buff[0].y * used->xsize + buff[0].x ] = label;
                *size_buff = 1; 
                break;
              }
          if( *size_buff == 1 ) break;
        } 
    }
  if( *size_buff == 0 ) return;  

  for( i=0; i < *size_buff; i++ )
    for( xx=buff[i].x-1; xx<=buff[i].x+1; xx++ )
      for( yy=buff[i].y-1; yy<=buff[i].y+1; yy++ )
        if( in_image( xx, yy, used->xsize, used->ysize ) )
          {
            adr = xx + yy * used->xsize;
            /* check if the point has already been used or visited */
            if( (used->data[adr] != USED) && (used->data[adr] != label) ) 
              {
                /* check if point inside the ellipse ring */                
                if( is_in_ell_ring( ering, xx, yy ) )
                  {
                    /* get angle normal to the ellipse */
                    theta = ellipse_normal_angle( (double)xx, (double)yy, foci);
                    /* check if gradient converges or diverges when computing 
                       the ideal angle of a point on the circle */
                    if( grad_dir == 0 )
                      {
		        if( theta > 0 ) theta = - ( M_PI - theta );
		        else theta = M_PI + theta;			
	              }
                    (*err) += norm_angle_diff( angles->data[adr], theta );
                    used->data[adr] = label;
                    buff[*size_buff].x = xx;
                    buff[*size_buff].y = yy;
                    ++(*size_buff); 
                  }
              }
          } 
  *pts = *size_buff;     
}


/*----------------------------------------------------------------------------*/
/** Iteration to improve NFA: count additive error and compute NFA. 
 */
static double ell_ring_nfa( PImageDouble angles, PImageInt used, Ring *ering,
                            double *foci, int x, int y, int grad_dir, 
                            Point *buff, int *size_buff, int label, 
                            double logNT_ell )
{
  int pts;
  double err;
  
  /* check parameters */
  if( (angles == NULL) || (angles->data == NULL) )
    error("ell_ring_nfa: invalid input 'angles' image.");
  if( (used == NULL) || (used->data == NULL) )
    error("ell_ring_nfa: invalid input 'used' image.");
  if( ering == NULL ) error("ell_ring_nfa: invalid input ring.");
  if( foci == NULL ) error("ell_ring_nfa: invalid input foci.");
  if( buff == NULL ) error("ell_ring_nfa: invalid local buffer.");
  if( size_buff == NULL ) 
    error("ell_ring_nfa: 'size_local_buff' must be non null.");
  if( label <= 1 ) error("ell_ring_nfa: forbidden label value.");

  /* count additive error */
  ell_count_error( angles, used, ering, foci, x, y, grad_dir, buff, size_buff,
                   label, &pts, &err );

  /* return computed NFA */
  return NFAc( pts, err, logNT_ell );
}


/*----------------------------------------------------------------------------*/
/** Try some ellipse ring variations to improve NFAc value.
 */
static double ell_ring_improve( PImageDouble angles, PImageInt used, 
                                Ring *ering, double *foci, int x, int y, 
                                int grad_dir, Point **local_buff, 
                                int *size_local_buff, Point **tmp_buff, 
                                int *label, double logNT_ell )
{ 
  double log_nfa, log_nfa_new;
  double delta = 0.5;
  int size_tmp_buff;

  /* check parameters */
  if( (angles == NULL) || (angles->data == NULL) )
    error("ell_ring_improve: invalid input 'angles' image.");
  if( (used == NULL) || (used->data == NULL) )
    error("ell_ring_improve: invalid input 'used' image.");
  if( ering == NULL ) error("ell_ring_improve: invalid input ring.");
  if( foci == NULL ) error("ell_ring_improve: invalid input foci.");
  if( local_buff == NULL ) error("ell_ring_improve: invalid local buffer.");
  if( tmp_buff == NULL ) error("ell_ring_improve: invalid tmp buffer.");
  if( size_local_buff == NULL ) 
    error("ell_ring_improve: 'size_local_buff' must be non null.");
  if( *label <= 1 ) error("ell_ring_improve: forbidden label value.");

  if ( ! check_ell_ring(ering) ) return -logNT_ell;
  (*label)++;
  *size_local_buff = 0;
  log_nfa = ell_ring_nfa( angles, used, ering, foci, x, y, grad_dir, 
                          *local_buff, size_local_buff, *label, logNT_ell );
     
  /* try to improve the ring, by reducing width; first reduce wmax */
  Ring r; 
  copy_ring( ering, &r );
  /*for( n=0; n<10; n++ )
    { */        
      while( r.wmax-delta > 0.0) 
        {
          r.wmax -= delta;
          size_tmp_buff = 0;
          (*label)++;
        
          log_nfa_new = ell_ring_nfa( angles, used, &r, foci, x, y, grad_dir,
                                      *tmp_buff, &size_tmp_buff, *label,
                                      logNT_ell );
          if( log_nfa_new > log_nfa )
            {
              ering->wmax -= delta;
              log_nfa = log_nfa_new;
              /* switch addresses of the pointers to the local and tmp buffer */
              /*swap_ptr_pts( local_buff, tmp_buff );
              (*size_local_buff) = size_tmp_buff;*/
            }            
        }
    /*} */

  /* now try to reduce wmin */
  copy_ring( ering, &r );
  /*for( n=0; n<10; n++ )
    { */         
      while( r.wmin + delta < 0.0) 
        {
          r.wmin += delta;
          size_tmp_buff = 0;
          (*label)++;
        
          log_nfa_new = ell_ring_nfa( angles, used, &r, foci, x, y, grad_dir,
                                      *tmp_buff, &size_tmp_buff, *label,
                                      logNT_ell );
          if( log_nfa_new > log_nfa )
            {
              ering->wmin += delta;
              log_nfa = log_nfa_new;
              /* switch addresses of the pointers to the local and tmp buffer */
              /*swap_ptr_pts( local_buff, tmp_buff );
              (*size_local_buff) = size_tmp_buff; */
            }            
        }  
   /* } */  
  if ( ! check_ell_ring(ering) ) return -logNT_ell;
  return log_nfa;
}   


/*----------------------------------------------------------------------------*/
/*------------------------------- BUILD OUTPUT -------------------------------*/
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/** A meaningful polygon has been found; store it in output. 
 */
void add_poly_out( Polygon **poly_out, PolyRect *poly, int poly_count, 
                   int *poly_count_max, int **labels, int lbl )
{
 
  if( poly_count>*poly_count_max )
    {
      *poly_out = (Polygon *) realloc ( *poly_out, (2*poly_count) *
                                                   sizeof(Polygon) );
   
      *labels = (int *) realloc ( *labels, (2*poly_count) * sizeof(int) );
      *poly_count_max = 2 * poly_count;
    }
  
  polyrect2polygon( poly, &( (*poly_out)[poly_count-1] ) );
  (*labels)[poly_count-1] = lbl;
}


/*----------------------------------------------------------------------------*/
/** A meaningful ellipse has been found; store it in output. 
 */
void add_ell_out( Ring **ell_out, Ring *ell, int ell_count, int *ell_count_max,
                  int **labels, int lbl )
{
  if( ell_count>*ell_count_max )
    {
      *ell_out = ( Ring *) realloc ( *ell_out, (2*ell_count) * sizeof(Ring) );
      *labels =  ( int *) realloc ( *labels, (2*ell_count) * sizeof(int) );
      *ell_count_max = 2 * ell_count;
    }
  
  copy_ring( ell, &( (*ell_out)[ell_count-1] ) );
  (*labels)[ell_count-1] = lbl;
}


/*----------------------------------------------------------------------------*/
/** Entry point in detector's code. Gets as input an image, and returns two 
    lists of geometric shapes: polygons and ellipses (circles are considered 
    particular ellipses, and segments are particular polygons). 
 */
void ELSDc( PImageDouble in, int *ell_count, Ring **ell_out, int **ell_labels,
            int *poly_count, Polygon **poly_out, int **poly_labels, 
            PImageInt out )
{
  double ang_th = 22.5;     /* gradient angle tolerance in degrees */
  double prec;              /* radian precision */
  double p;
  double rho = 5.0;         /* pixels with gradient magnitude under this 
                               threshold are not considered, as their (weak) 
                               orientation might be too affected by noise  */
  int grad_dirc;            /* direction of the gradient wrt center of the 
                               circle */
  int grad_dire;            /* direction of the gradient wrt center of the 
                               ellipse */
  double density_th = 0.7;  /* minimum density */
  void *mem_p;              /* parameters used for sorting */ 
  int n_bins = 1024;        /* (binning) pixels according to */  
  double max_grad = 255.0;  /* to their gradient magnitude */
  CoordList *list_p;        /* list of sorted pixels */
  unsigned int xsize,ysize; /* image size */
  double mlog10eps = 0.0;   /* minus log of the number of accepted false 
                               positives ( eps=1 ) */
  int min_size_seg;         /* minimum size of a region that could contain a 
                               meaningful segment */
  int min_size_ell;         /* minimum size of a region that could contain a 
                               meaningful ellipse */
  int min_size_poly;        /* minimum size of a region that could contain a 
                               meaningful ellipse */
  double msize;             /* max value allowed for circle/ellipse radius */
  double logNT_ell;         /* log of the number of tests for elliptical arcs */
  double logNT_seg;         /* log of the number of tests for line segments */ 
  double logNT_poly;        /* log of the number of tests for polygons */
  double log_nfa, best_nfa; /* NFA values using the continuous */
  void *best_feature;       /* pointer to best feature found so far */
  char best_type = '0';     /* type of best feature; 'p' for polygon, 
                               'e' for ellipse */
  PImageDouble imgauss;     /* smooth version of the original image, used during
                               candidate selection */
  PImageDouble angles;      /* gradient angles smooth image */
  PImageDouble angles0;     /* gradient angles original image */
  PImageDouble gradmag;     /* gradient magnitude */
  PImageDouble gradx, grady;/* gradient orientations on Ox and Oy */
  PImageInt used;           /* image to mark already used points */
  PolyRect *poly;           /* parameters of polygon as list of rectangles */
  PolyRect *seg;            /* one rectangle is a particular polygon */
  Ring ering, cring;        /* ring parameters */
  double *buff_fit;         /* buffer zone needed for fitting */
  int size_buff_fit;        /* size of the currently allocated memory for 
                               the buffer */
  Point *reg;               /* list of points gathered during curve growing */
  int reg_size;
  Point *best_buff;         /* buffer to store the points of the most meaningful
                               primitive found so far (among segment, polygon, 
                               ellipse, circle) in one iteration */
  int size_best_buff;
  Point *new_buff;          /* buffer to store the points of a new primitive; 
                               if the new primitive is more meaningful than 
                               the most meaningful primitive found so far, 
                               deem new primitive as most meaningful by 
                               switching the addresses with best_buff   */
  int size_new_buff;        
  Point *tmp_buff;          /* buffer to store the points of a new version of a
                               primitive, computed during improve operation; 
                               if the new version is more meaningful than 
                               the most meaningful version found so far, 
                               deem new version as most meaningful by 
                               switching the addresses with new_buff   */
  int ell_count_max = 0;
  int poly_count_max = 0;
  int adr;
  double cparam[5];
  double eparam[5]; 
  double pext1[2], pext2[2], ang0;
  int label = 2;            /* 0 is NOTUSED, 1 is USED, so 'label' must go 
                               upward from 2 */              
  int label_meaningful = 0;
  double spir;
  double foci[4];
  /*int i;*/

  /* check parameters */
  if( (in == NULL) || (in->data == NULL) ) error("ELSDc: invalid input image.");
  if( (out == NULL) || (out->data == NULL) ) 
    error("ELSDc: invalid output label image.");

  /* allocate and init variables */
  prec = M_PI*ang_th/180.0;
  p = ang_th/180.0;
  xsize = in->xsize;
  ysize = in->ysize;
  
  /* allocate matrix to store elements of the design matrix used in conic
     fitting */
  buff_fit = (double*) malloc ( sizeof (double) );
  size_buff_fit = 1;

  /* allocate buffers to store point coordinates */
  reg = (Point *) malloc( xsize * ysize * sizeof(Point) );
  if( reg == NULL ) error("ELSDc: not enough memory.");
  best_buff = (Point *) malloc( xsize * ysize * sizeof(Point) );
  if( best_buff == NULL ) error("ELSDc: not enough memory.");
  new_buff = (Point *) malloc( xsize * ysize * sizeof(Point) );
  if( new_buff == NULL ) error("ELSDc: not enough memory.");
  tmp_buff = (Point *) malloc( xsize * ysize * sizeof(Point) );
  if( tmp_buff == NULL ) error("ELSDc: not enough memory.");

  used = new_PImageInt_ini( xsize, ysize, NOTUSED );

  /* Number of tests for ellipse arcs. We must count (3*xsize)*(3*ysize) 
     possible values for the center, as center can be outside the image too,
     (xsize*ysize) possible values for the two axes, sqrt(xsize*ysize) possible
     values for orientation, sqrt(xsize*ysize) possible values for width, 
     (xsize*ysize) possible values for delimiting angles. This makes 
     9*(xsize*ysize)^4. Since we have polygons and ellipses, i.e. two families,
     we must multiply this by 2. So NT_ell = 2*9*(xsize*ysize)^4. 
     In log, we get:   */
  logNT_ell = 4.0 * ( log10( (double)xsize ) + log10( (double)ysize ) ) + 
                      LOG10_9 + LOG10_2; 
  min_size_ell = (int)( (-logNT_ell-mlog10eps) / log10(p) );

  logNT_seg = 2.5 * (log10((double)xsize) + log10((double)ysize) ) 
                     + LOG10_2;
  min_size_seg = (int)( (-logNT_seg - mlog10eps) / log10(p) );

  /* perform gaussian smoothing */
  imgauss = gaussian_sampler( in, 1.0, 0.6 );

  /* compute gradient magnitude and orientation  for the smooth image */
  img_gradient_sort( imgauss, rho, &list_p, &mem_p, n_bins, max_grad, &angles, 
                     &gradmag, &gradx, &grady );

  /* compute gradient orientation for the original image */
  angles0 = img_gradient_angle( in, rho );
 
  /* input image not needed any more; free it */
  free_PImageDouble(in);
  
  msize = max(xsize,ysize);

  poly = new_polyrect();
  seg = new_polyrect();

  
  /* begin primitive detection */
  for( ; list_p; list_p = list_p->next )
    { 
      adr = list_p->y * used->xsize + list_p->x;
      /* only pixels that haven't been used already can be seed pixels */
      if( used->data[adr] == USED ) continue;
      
      /* only pixels with defined gradient can be seed pixels */
      if( angles->data[adr] == NOTDEF ) continue;
      
      /* only pixels that are local maxima in the image can be seed pixels */  
      /*if( !local_max( imgauss, list_p->x, list_p->y ) ) continue; */
     
      /* init intra-iteration variables */ 	
      log_nfa = best_nfa = mlog10eps;
      best_feature = NULL; 
      best_type = '0';
      size_best_buff = size_new_buff = 0;
      label++;
      spir = 0.0; 
      clear_polyrect(poly);
      clear_polyrect(seg);
          
      /* set initial seed point */
      reg[0].x = list_p->x; 
      reg[0].y = list_p->y;
      reg_size = 1;

      /* Gather points aligned along a convex polygon */
      if( !curve_grow( gradmag, angles, used, reg, &reg_size, density_th,
                       prec, poly, &label, pext2, pext1, &spir ) ) 
        continue;

      
      /* Validate polygon */
      /* Compute number of tests for polygon; it must be done for each 
         polygon, as it depends on the number of segments in the 
         polygon */
      logNT_poly = ( poly->dim + 1.5 ) * ( log10((double)xsize ) + 
                     log10( (double)ysize) ) + ( poly->dim+1.0 ) * LOG10_2;
      /* polygon's best validation score */
      min_size_poly =(int)((-logNT_poly - mlog10eps)/log10(p));
      
      if( reg_size > min_size_poly )
        {
          best_nfa = poly_improve( angles0, used, poly, &best_buff, 
                                   &size_best_buff, &tmp_buff, logNT_poly );
          best_feature = (void*)poly;
          best_type = 'p';
        }         

      /* Estimate a single segment on the entire set and validate; 
         this prevents from considering very small regions at the ends 
         of a large segment as new segments. */
      get_seg( angles, gradmag, used, seg, reg, reg_size, prec );
      log_nfa = poly_improve( angles0, used, seg, &new_buff, 
                              &size_new_buff, &tmp_buff, logNT_seg );
    
      /* if more meaningful than previous, update best feature */
      if( log_nfa > best_nfa && size_new_buff>min_size_seg ) 
        {
          best_nfa = log_nfa;
          best_feature = (void*)seg;
          best_type = 'p';
          swap_ptr_pts( &best_buff, &new_buff );
          size_best_buff = size_new_buff;
        }

      /* Estimate circle and ellipse parameters (centre, axes, orientation)
         on the gathered pixels */
      if( reg_size > min_size_ell ) 
        conic_fit( gradx, grady, reg, reg_size, &buff_fit, &size_buff_fit, 
                   cparam, eparam );

      /* Get parameters of the circle ring and compute validation score */
      ang0 = angles->data[ reg[0].y*angles->xsize + reg[0].x ];
      if( get_ring( reg, reg_size, ang0, cparam, pext1, pext2, CIRCLE,
                    spir, &cring, &grad_dirc, foci, msize ) )
        {
          log_nfa = circ_ring_improve( angles0, used, &cring, reg[0].x, 
                                       reg[0].y, grad_dirc, &new_buff, 
                                       &size_new_buff, &tmp_buff, &label,
                                       logNT_ell );
          /* if more meaningful than previous, update best feature */
          if( log_nfa > best_nfa && size_new_buff>min_size_ell ) 
            {
              best_nfa = log_nfa;
              best_feature = (void*)(&cring);
              best_type = 'c';
              swap_ptr_pts( &best_buff, &new_buff );
              size_best_buff = size_new_buff;
            }
        } /* ! get cring */
   
      /* Get parameters of the ellipse ring and compute validation score */
      if( get_ring( reg, reg_size, ang0, eparam, pext1, pext2, ELLIPSE,
                    spir, &ering, &grad_dire, foci, msize ) )
        {
          log_nfa = ell_ring_improve( angles0, used, &ering, foci, reg[0].x, 
                                      reg[0].y, grad_dire, &new_buff, 
                                      &size_new_buff, &tmp_buff, &label,
                                      logNT_ell ); 
          /* if more meaningful than previous, update best feature */
          if( log_nfa > best_nfa && size_new_buff>min_size_ell ) 
            {
              best_nfa = log_nfa;
              best_feature = (void*)(&ering);
              best_type = 'c';
              swap_ptr_pts( &best_buff, &new_buff );
              size_best_buff = size_new_buff;
            }
        }  /* ! get ering */

      /* At this point, we have the most meaninful feature at this iteration. 
         If it is meaningful compared to eps, then store */   
      if( best_nfa <= mlog10eps ) continue;
      label_meaningful++;
      if( best_type == 'c' )
        {
          (*ell_count)++;
          add_ell_out( ell_out, (Ring*)best_feature, *ell_count, &ell_count_max, 
                       ell_labels, label_meaningful );
          /*Ring *rr = (Ring*)best_feature;*/
          /* free(rr); */
        }
      else if( best_type == 'p' ) 
        {
          (*poly_count)++;
          add_poly_out( poly_out, (PolyRect*)best_feature, *poly_count, 
                        &poly_count_max, poly_labels, label_meaningful );
        }
 

      /* Mark as USED the pixels of the best feature in 'used' */
      mark_img_pts( used, best_buff, 0, size_best_buff, USED );

      /* Mark label of the best feature in output label image */
      mark_img_pts( out, best_buff, 0, size_best_buff, label_meaningful );  
   
    }  /* ! for loop */

  /* free useless structures and arrays */		    
  free_PImageDouble(imgauss);
  free_PImageDouble(gradx); 
  free_PImageDouble(grady);
  free_PImageDouble(gradmag); 
  free_PImageDouble(angles);
  free_PImageDouble(angles0);  
  free_PImageInt(used);
  free(reg); 
  free(best_buff);
  free(new_buff);
  free(tmp_buff);
  free(buff_fit);  
  free(poly->rectlist); free(poly); 
  free(seg->rectlist); free(seg); 
  free(mem_p);
}
