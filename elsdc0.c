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
#include "elsdc.h"


/*----------------------------------------------------------------------------*/
/** Continuous formulation of the number of false alarms 
    NFAC = NT * k^n / n!
    log(n!) is bounded by Stirling's approximation:
       n! >= sqrt(2pi) * n^(n+0.5) * exp(-n)
     then, log10(NFA) <= log10(NT) + n*log10(k) - log10(latter expansion) 
 */ 
#define 1_2_LOG10_M_2__PI  0.39908993418
#define LOG10_EXP1         0.43429448190
static double NFAc( int n, double k, double logNT )
{
  /* check parameters */
  if( n<= 0) error("NFAc: 'n' must be strictly positive.");

  double logNFAC = 0.0;
  if( k>0 ) 
    logNFAC = logNT + n * log10(k) - 1_2_LOG10_M_2__PI - 
              (n+0.5) * log10(n) + n * LOG10_EXP1 ;

  return -logNFAC;
}


/*----------------------------------------------------------------------------*/
/** Determine if the gradient converges or diverges to/from the centre.
    Return TRUE if gradient diverges, FALSE if it converges.
    Circle case.
 */
static int convergent_gradient_circ( int px, int py, double xc, double yc )
{
  double a;
  int dir = FALSE;
  a = atan2( py-yc, px-xc );
  if( angle_diff( a, angles->data[ py * angles->xsize + px ] ) < M_1_2_PI )
    dir = TRUE;
  return dir;
}


/*----------------------------------------------------------------------------*/
/** Determine if the gradient converges or diverges to/from the centre.
    Return TRUE if gradient diverges, FALSE if it converges.
    Ellipse case.
 */
static int convergent_gradient_ell( int px, int py, double *foci )
{
  double a;
  int dir = FALSE;

  /* check parameters */
  if( foci == NULL) error("int_ext_e: invalid ellipse.");
  a = ellipse_normal_angle( (double)px, (double)py, foci );
  if( angle_diff( a, angles->data[ py * angles->xsize + px ] ) < M_1_2_PI )
    dir = TRUE;
  return dir;
}


/*----------------------------------------------------------------------------*/
/** Group in 'reg' neighbour pixels sharing the same orientation up to 
    precision 'prec'.
    reg_index: 0.....start.....idx_buf (end)
    (start...idx_buf) is the buffer zone (>=1) of the list; the orientations of 
    the pixels in this zone are used to initialize reg_angle. 
 */
static void region_grow( PImageDouble angles, PImageInt used, Point *reg, 
                         int start, int idx_buff, int *end, int label, 
                         double prec, double *reg_angle )
{
  double sumdx = 0.0, sumdy = 0.0;
  double ang;
  int xx, yy, i;
  int adr;

  /* check parameters */
  if( (angles == NULL) || (angles->data == NULL) )
    error("region_grow: invalid image 'angles'.");
  if( (used == NULL) || (used->data == NULL) )
    error("region_grow: invalid image 'used'.");
  if( reg == NULL ) error("region_grow: invalid input region.");
  if( start >= *end ) error("region_grow: invalid indexes.");
  if( prec < 0.0 ) error("region_grow: 'prec' must be positive.");
  
  /* Initialise reg_angle with the value of the pixels 
     orientations already added to the region. First time when it is called,  */
  for( i=start; i<idx_buff; i++ )
    {
      xx = reg[i].x; yy = reg[i].y;
      adr = xx + yy*angles->xsize;
      sumdx += cos( angles->data[adr] );
      sumdy += sin( angles->data[adr] );
      *reg_angle = atan2( sumdy, sumdx );
      used->data[adr] = label;
    }
  
  /* try neighbors as new region points */
  for( i=start; i<*end; i++ )
    for( xx=reg[i].x-1; xx<=reg[i].x+1; xx++ )
      for( yy=reg[i].y-1; yy<=reg[i].y+1; yy++ )
        if( in_image( xx, yy, used ) )
          {
            adr = yy * angles->xsize + xx;
            if( (used->data[adr] != USED) && (used->data[adr] != label) )
              {
                ang = angles->data[adr];
                if( (ang != NOTDEF ) && (is_similar( ang, *reg_angle, prec )) )
                  {
                    /* add point and mark as already visited in this call */
                    used->data[adr] = label;            
                    reg[*end].x = xx;
                    reg[*end].y = yy;
                    ++(*end);             
                    /* update region's angle */
                    sumdx += cos(ang);
                    sumdy += sin(ang);
                    *reg_angle = atan2( sumdy, sumdx );
                  }
         
              }
          }
}


/*----------------------------------------------------------------------------*/
/*--------- FUNCTIONS FOR RECTANGLE ESTIMATION AND VALIDATION ----------------*/
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/** Estimate the angle of the rectangle covering the pixels listed in 'reg' 
    between 'start' and 'end' indexes. 'cx' and 'cy' are the centroid 
    coordinates of the point set. The modulus of the gradient magnitude 
    contained in 'gradmag' is used to weight pixels' contributions. 
 */
static double get_theta( PDoubleImage gradmag, Point *reg, int start, int end, 
                         double cx, double cx, double reg_angle, double prec )
{
  double lambda1, lambda2, tmp, theta, weight, sum;
  double Ixx = 0.0;
  double Iyy = 0.0;
  double Ixy = 0.0;
  int i;
  
  /* check parameters */
  if( reg == NULL ) error("get_theta: invalid input region.");
  if( (gradmag == NULL) || (gradmag->data == NULL) )
    error("get_theta: invalid image 'gradmag'.");
  if( start >= end ) error("get_theta: invalid indexes.");
  if( prec < 0.0 ) error("get_theta: 'prec' must be positive.");

  /*----------- theta ---------------------------------------------------*/
  /*
      Region inertia matrix A:
         Ixx Ixy
         Ixy Iyy
      where
        Ixx = \sum_i y_i^2
        Iyy = \sum_i x_i^2
        Ixy = -\sum_i x_i y_i

      lambda1 and lambda2 are the eigenvalues, with lambda1 >= lambda2.
      They are found by solving the characteristic polynomial
      det(\lambda I - A) = 0.

      To get the line segment direction we want to get the eigenvector of
      the smaller eigenvalue. We have to solve a,b in:
        a.Ixx + b.Ixy = a.lambda2
        a.Ixy + b.Iyy = b.lambda2
      We want the angle theta = atan(b/a). I can be computed with
      any of the two equations:
        theta = atan( (lambda2-Ixx) / Ixy )
      or
        theta = atan( Ixy / (lambda2-Iyy) )

      When |Ixx| > |Iyy| we use the first, otherwise the second
      (just to get better numeric precision).
   */
  sum = 0.0;
  for( i=start; i<end; i++ )
    {
      weight = gradmag->data[ reg[i].x + reg[i].y * gradmag->xsize ];
      Ixx += ( (double) reg[i].y - cy ) * ( (double) reg[i].y - cy ) * weight;
      Iyy += ( (double) reg[i].x - cx ) * ( (double) reg[i].x - cx ) * weight;
      Ixy -= ( (double) reg[i].x - cx ) * ( (double) reg[i].y - cy ) * weight;
      sum += weight;
    }
  if( sum <= 0.0 ) error("get_theta: weights sum less or equal to zero.");
  Ixx /= sum;
  Iyy /= sum;
  Ixy /= sum;
  lambda1 = ( Ixx + Iyy + sqrt( (Ixx-Iyy)*(Ixx-Iyy) + 4.0*Ixy*Ixy ) ) / 2.0;
  lambda2 = ( Ixx + Iyy - sqrt( (Ixx-Iyy)*(Ixx-Iyy) + 4.0*Ixy*Ixy ) ) / 2.0;
  if( fabs(lambda1) < fabs(lambda2) )
    {
      fprintf(stderr,"Ixx %g Iyy %g Ixy %g lamb1 %g lamb2 %g - lamb1 < lamb2\n",
              Ixx, Iyy, Ixy, lambda1, lambda2 );
      tmp = lambda1;
      lambda1 = lambda2;
      lambda2 = tmp;
    }

  if( fabs(Ixx) > fabs(Iyy) )
    theta = atan2( lambda2-Ixx, Ixy );
  else
    theta = atan2( Ixy, lambda2-Iyy );

  /* The previous procedure doesn't care about orientation,
     so it could be wrong by 180 degrees. Here is corrected if necessary. */
  if( angle_diff( theta-M_1_2_PI, reg_angle ) > prec ) theta += M_PI;

  /* If the estimated orientation is still very different wrt reg_angle
     (region is nearly square), return the orthogonal to original reg_angle. */
  if( angle_diff( theta-M_1_2_PI, reg_angle ) > prec ) 
    theta = reg_angle + M_1_2_PI;

  while( theta >  M_PI ) theta -= M_2__PI;
  while( theta < -M_PI ) theta += M_2__PI;

  return theta;
}


/*----------------------------------------------------------------------------*/
/** Compute a rectangle that covers a connected region of points. A rectangle is
    defined by end points, center, orientation, width. Additional fields 
    (precision of approximation) are useful for validation. 
 */
static void region2rect( PDoubleImage gradmag, Point *reg, int start, int end,
                         double reg_angle, Rectangle *rec, double prec )
{
  double cx, cy, dx, dy, l, w, theta;
  double weight,sum,l_min,l_max,w_min,w_max;
  int i, adr;

  /* check parameters */
  if( (gradmag == NULL) || (gradmag->data == NULL) )
    error("region2rect: invalid image 'gradmag'.");
  if( reg == NULL ) error("region2rect: invalid input region.");
  if( start >= end ) error("region2rect: invalid indexes."); 
  if( rec == NULL ) error("region2rect: invalid 'rec'.");

  /* center */
  cx = cy = sum = 0.0;
  for( i=start; i<end; i++ )
    {
      adr = reg[i].x + reg[i].y * gradmag->xsize;
      weight = gradmag->data[adr];
      cx += (double) reg[i].x * weight;
      cy += (double) reg[i].y * weight;
      sum += weight;
    }

  cx /= sum;
  cy /= sum;

  /* orientation */  
  theta = get_theta( gradmag, reg, start, end, cx, cy, prec );

  /* length and width */
  dx = cos(theta);
  dy = sin(theta);
  l_min = l_max = w_min = w_max = 0.0;
  for( i=start; i<end; i++ )
    {
      l =  ( (double) reg[i].x - cx) * dx + ( (double) reg[i].y - cy) * dy;
      w = -( (double) reg[i].x - cx) * dy + ( (double) reg[i].y - cy) * dx;

      if( l > l_max ) l_max = l;
      if( l < l_min ) l_min = l;
      if( w > w_max ) w_max = w;
      if( w < w_min ) w_min = w;
    }

  /* store values */
  rec->x1 = cx + l_min * dx;
  rec->y1 = cy + l_min * dy;
  rec->x2 = cx + l_max * dx;
  rec->y2 = cy + l_max * dy;
  rec->width = w_max - w_min;
  rec->wmin = w_min;
  rec->wmax = w_max;
  rec->xc = cx;
  rec->yc = cy;
  rec->theta = theta;
  rec->dx = dx;
  rec->dy = dy;
  rec->prec = prec;

  /* correct if rectangle too thin */
  if( rec->width < 1.0 )
    { 
      rec->width = 1.0; rec->wmin = -0.5; rec->wmax = 0.5;
    }
}


/*----------------------------------------------------------------------------*/
/** Reduce the region size, by eliminating the points far from the starting 
    point, until that leads to rectangle with the right density of region points
    or to discard the region if too small.
 */
static int reduce_region_radius( PImageDouble gradmag, PImageDouble angles,
                                 PImageInt used, Point *reg, int start, 
                                 int *end, double reg_angle, Rectangle *rec, 
                                 int x, int y, double density_th, double prec )
{
  double density, rad1, rad2, rad;
  int i;

  /* check parameters */
  if( (gradmag == NULL) || (gradmag->data == NULL) )
    error("reduce_region_radius: invalid image 'gradmag'.");
  if( (angles == NULL) || (angles->data == NULL) )
    error("reduce_region_radius: invalid image 'angles'.");
  if( (used == NULL) || (used->data == NULL) )
    error("reduce_region_radius: invalid image 'used'.");
  if( reg == NULL ) error("reduce_region_radius: invalid input region.");
  if( start >= *end ) error("reduce_region_radius: invalid indexes."); 
  if( rec == NULL ) error("reduce_region_radius: invalid 'rec'.");
  if( prec < 0.0 ) error("reduce_region_radius: 'prec' must be positive.");

  /* compute region points density */
  density = (double) (*end - start) /
                     (dist(rec->x1,rec->y1,rec->x2,rec->y2) * rec->width);

  if( density >= density_th ) return TRUE;

  /* compute region radius */
  rad1 = dist( (double)x, (double)y, rec->x1, rec->y1 );
  rad2 = dist( (double)x, (double)y, rec->x2, rec->y2 );
  rad = rad1 > rad2 ? rad1 : rad2;

  while( density < density_th )
    {
      rad *= 0.75;

      /* remove points from the region and update 'used' map */
      for( i=start; i<*end; i++ )
        if( dist( x, y, (double) reg[i].x, (double) reg[i].y ) > rad )
          {
            /* point not kept, mark it as NOTUSED */
            used->data[ reg[i].x + reg[i].y * used->xsize ] = NOTUSED;
            /* remove point from the region */
            reg[i].x = reg[*end-1].x; /* if i==*end-1 copy itself */
            reg[i].y = reg[*end-1].y;
            --(*end);
            --i; /* to avoid skipping one point */
          }
      /* reject if the region is too small.
         2 is the minimal region size for 'region2rect' to work. */

      if( (*end - start) < 2 ) return FALSE;

      /* re-compute rectangle */
      region2rect( gradmag, reg, start, *end, reg_angle, rec, prec );

      /* re-compute region points density */
      density = (double) ( *end - start) /
                ( dist( rec->x1, rec->y1, rec->x2, rec->y2 ) * rec->width );     
    }
  return TRUE;
}


/*----------------------------------------------------------------------------*/
/** Refine a rectangle. For that, an estimation of the angle tolerance is
    performed by the standard deviation of the angle at points near the
    region's starting point. Then, a new region is grown starting from the
    same point, but using the estimated angle tolerance.
    If this fails to produce a rectangle with the right density of
    region points, 'reduce_region_radius' is called to try to
    satisfy this condition.
 */
static int refine( PImageDouble gradmag, PImageDouble angles, PImageInt used, 
                   Point *reg, int start, int *end, int idx_buff, int *label,
                   Rectangle *rec, double density_th, double prec )
{
  double angle, ang_d, mean_angle, tau, density, ang_c, sum, s_sum, reg_angle;
  double gradmax;
  int i, n, xx, yy;
  int xc, yc;

  /* check parameters */
  if( (gradmag == NULL) || (gradmag->data == NULL) ) 
    error("refine: invalid image 'gradmag'.");
  if( (angles == NULL) || (angles->data == NULL) ) 
    error("refine: invalid image 'angles'.");
  if( (used == NULL) || (used->data == NULL) ) 
    error("refine: invalid image 'used'.");
  if( reg == NULL ) error("refine: invalid input region.");
  if( start >= *end ) error("refine: invalid indexes."); 
  if( *label <= 1 ) error("refine: forbidden label value.");
  if( prec < 0.0 ) error("refine: 'prec' must be positive.");
  if( rec == NULL ) error("refine: invalid pointer 'rec'.");
 
  /* compute region points density */
  density = (double) (*end - start) /
            ( dist(rec->x1,rec->y1,rec->x2,rec->y2) * rec->width );
 
  if( density >= density_th ) return TRUE;
  
  /*------ First try: reduce angle tolerance ------*/

  /* compute the new mean angle and tolerance */
  gradmax = 0.0;
  xc = reg[start].x;
  yc = reg[start].y;

  for( i=start; i<idx_buff; i++ )
    {
      xx = reg[i].x; yy = reg[i].y;
      if( gradmag->data[yy*gradmag->xsize+xx] > gradmax )
        {
          xc = xx;
          yc = yy;
          gradmax = gradmag->data[yy*gradmag->xsize+xx];
        }
    }
  ang_c = angles->data[ xc + yc * angles->xsize ];
  sum = s_sum = 0.0;
  n = 0;
  for( i=start; i<*end; i++ )
    {
      used->data[ reg[i].x + reg[i].y * used->xsize ] = NOTUSED;
      if( dist( (double)xc,(double)yc,(double)reg[i].x,(double)reg[i].y ) 
          < rec->width )
        {
          angle = angles->data[ reg[i].x + reg[i].y * angles->xsize ];
          ang_d = angle_diff_signed(angle,ang_c);
          sum += ang_d;
          s_sum += ang_d *ang_d;
          ++n;
        }
    }
  mean_angle = sum / (double) n;
  /* tau = 2 * standard deviation */
  tau = max(2.0 * sqrt( (s_sum - 2.0 * mean_angle * sum) / (double) n
                         + mean_angle * mean_angle ),0.2); 
  tau = min(tau,prec);
  /*tau = 2.0 * sqrt( (s_sum - 2.0 * mean_angle * sum) / (double) n
                      + mean_angle*mean_angle );*/ /* 2 * standard deviation */

 
  /* find a new region from the same starting point and new angle tolerance */
  *end = idx_buff;
  (*label)++;
  region_grow(angles, used, reg, start, idx_buff, end, *label, tau, &reg_angle);

  /* if the region is too small, reject */
  if( *end - start <= 2 ) return FALSE;

  /* re-compute rectangle */
  region2rect( gradmag, reg, start, *end, reg_angle, rec, prec );

  /* re-compute region points density */
  density = (double) (*end - start) /
                     ( dist(rec->x1,rec->y1,rec->x2,rec->y2) * rec->width);

  /*------ Second try: reduce region radius ------*/
  if( density < density_th )
    return reduce_region_radius( gradmag, angles, used, reg, start, end, 
                                 reg_angle, rec, xc, yc, density_th, prec );
  return TRUE;
}


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

  for( i=ini_RectIter(rec); !end_RectIter(i); inc_RectIter(i) )
    if( in_image( i->x, i->y, used ) ) 
      {
        ++(*pts);
        adr = (i->y) * used->xsize + (i->x);
        if( (used->data[adr] == USED) || (angles->data[adr] == NOTDEF) )
          (*alg) += 1.0;
        else
          (*alg) += norm_angle_diff( angles->data[adr], theta );

        buff[*size_buff].x = i->x;
        buff[*size_buff].y = i->y;
        ++(*size_buff);
      }

  free_RectIter(i);
}


/*----------------------------------------------------------------------------*/
/** Compute a rectangle's NFA value.
 */
static double rect_nfa( PImageDouble angles, PImageInt used, Rectangle *rec,
                        Point *buff, int *size_buff, double logNT )
{
  int pts; 
  double err;

  /* check parameters */
  if( (angles == NULL) || (angles->data == NULL) )
    error("rect_count_error: invalid input 'angles' image.");
  if( (used == NULL) || (used->data == NULL) )
    error("rect_count_error: invalid input 'used' image.");
  if( rec == NULL ) error("rect_count_error: invalid input rectangle.");
  if( buff == NULL ) error("rect_count_error: invalid buffer.");
  if( size_buff == NULL ) 
    error("rect_count_error: size_buff must be non null.");
  
  /* count additive error */
  rect_count_error( angles, used, rec, &pts, &err, buff, size_buff ); 

  /* return computed NFA */
  return NFAc( pts, alg, logNT );
}


/*----------------------------------------------------------------------------*/
/** Try some rectangle variations to improve NFAc value.
 */
static double rect_improve( PImageDouble angles, PImageInt used, Rectangle *rec, 
                            Point **local_buff, int *size_local_buff, 
                            Point **tmp_buff, double logNT )
{
  Rectangle r;
  double log_nfa,log_nfa_new;
  double delta = 0.5;
  double delta_2 = delta / 2.0;
  int n;
  int size_tmp_buff;
  Point *tmp;

  /* check parameters */
  if( (angles == NULL) || (angles->data == NULL) )
    error("rect_improve: invalid input 'angles' image.");
  if( (used == NULL) || (used->data == NULL) )
    error("rect_improve: invalid input 'used' image.");
  if( rec == NULL ) error("rect_improve: invalid input rectangle.");
  if( *local_buff == NULL ) error("rect_improve: invalid local buffer.");
  if( size_local_buff == NULL ) 
    error("rect_improve: 'size_local_buff' must be non null.");
  if( *tmp_buff == NULL ) error("rect_improve: invalid tmp buffer.");
 
  *size_local_buff = 0;
  log_nfa = rect_nfa( angles, used, rec, *local_buff, size_local_buff, logNT );

  /* try to reduce dmax */
  rect_copy( rec, &r );     
  for( n=0; n<10; n++ )
    {
      if( r.wmax - delta > 0.0 )
        { 
          r.wmax -= delta;
          size_tmp_buff = 0;
          log_nfa_new = rect_nfa( angles, used, &r, *tmp_buff, &size_tmp_buff,
                                  logNT );
          if( log_nfa_new > log_nfa )
            {
              rect_copy(&r,rec);
              log_nfa = log_nfa_new;
              /* switch addresses of the pointers to the local and tmp buffer */
              swap_ptr_pts( local_buff, tmp_buff );
              (*size_local_buff) = size_tmp_buff;
            }
        }
    }

  /* now try to reduce dmin */
  rect_copy( rec, &r );
  for( n=0; n<10; n++ )
    { 
      if( r.wmin + delta < 0.0 )
        {
          r.wmin += delta;
          size_tmp_buff = 0;
          log_nfa_new = rect_nfa( angles, used, &r, *tmp_buff, &size_tmp_buff,
                                  logNT );
          if( log_nfa_new > log_nfa )
            {
              rect_copy( &r, rec );
              log_nfa = log_nfa_new;
              /* switch addresses of the pointers to the local and tmp buffer */
              swap_ptr_pts( local_buff, tmp_buff );
              (*size_local_buff) = size_tmp_buff;
            }
        }  
    }  
  return log_nfa;
}


/*----------------------------------------------------------------------------*/
/** Given a list of points, estimate rectangle and compute best validation 
    score.
 */
static double valid_seg( PImageDouble angles, PImageDouble gradmag, 
                         PImageInt used, Rectangle *rec, Point *reg, 
                         int reg_size, Point **local_buff, int *size_local_buff, 
                         Point **tmp_buff, double prec, double logNT_seg )
{
  double sumdx = 0.0, sumdy = 0.0;
  int i;
  double nfa_seg;
  int adr;

  /* check parameters */
  if( (angles == NULL) || (angles->data == NULL) )
    error("valid_seg: invalid input 'angles' image.");
  if( (gradmag == NULL) || (gradmag->data == NULL) )
    error("valid_seg: invalid input 'angles' image.");
  if( (used == NULL) || (used->data == NULL) )
    error("valid_seg: invalid input 'used' image.");
  if( rec == NULL ) error("valid_seg: invalid input rectangle.");
  if( reg == NULL ) error("valid_seg: invalid input region.");
  if( reg_size <= 0 ) error("valid_seg: 'reg_size' must be strictly positive.");
  if( *local_buff == NULL ) error("valid_seg: invalid local buffer.");
  if( buff_sz == NULL ) 
    error("valid_seg: local buffer size 'sz' must be non null.");
  if( *tmp_buff == NULL ) error("valid_seg: invalid tmp buffer.");
  
  /* compute the mean angle of the whole region */
  for( i=0; i<reg_size; i++ )
    {
      adr = reg[i].x+reg[i].y*angles->xsize;
      sumdx += cos( angles->data[adr] );
      sumdy += sin( angles->data[adr] );
    }
  reg_angle = atan2( sumdy, sumdx );

  /* estimate rectangle */
  region2rect( gradmag, reg, 0, reg_size, reg_angle, rec, prec ); 

  /* compute best validation score */
  nfa_seg = rect_improve( angles, used, rec, local_buff, size_local_buff, 
                          tmp_buff, logNT_seg ); 
  return nfa_seg;     
}


/*----------------------------------------------------------------------------*/
/*--------- FUNCTIONS FOR POLYGON ESTIMATION AND VALIDATION ------------------*/
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/** Compute additive normalised per-pixel error inside a polygon, and stack
    visited points so they can be marked as 'used' if the polygon is the most 
    meaningful primitive.
 */
static void poly_count_error( PImageDouble angles, PImageInt used, 
                              PolyRect *poly, int *pts, double *err, 
                              Point *buff, int *size_buff ) 
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
      poly->rectlist[i].wmin = poly->wmin;
      poly->rectlist[i].wmax = poly->wmax;
      rect_count_error( angles, used, &(poly->rectlist[i]), &pts_rect, 
                        &err_rect, &(buff[size_buff]), size_buff );
      (*pts) += pts_rect;
      (*err) += err_rect;
    }  
}


/*----------------------------------------------------------------------------*/
/** Iteration to improve NFA: count additive error and compute NFA. 
 */
static double poly_nfa( PImageDouble angles, PImageInt used, PolyRect *poly, 
                        Point *buff, int *size_buff, double logNT_poly )
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
  poly_count_error( angles, used, poly, &pts, &err, buff, size_buff );
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
  int n;
  int size_tmp_buff;
  Point *tmp;

  /* check parameters */
  if( (angles == NULL) || (angles->data == NULL) )
    error("poly_improve: invalid input 'angles' image.");
  if( (used == NULL) || (used->data == NULL) )
    error("poly_improve: invalid input 'used' image.");
  if( rec == NULL ) error("poly_improve: invalid input rectangle.");
  if( *local_buff == NULL ) error("poly_improve: invalid local buffer.");
  if( size_local_buff == NULL ) 
    error("poly_improve: 'size_local_buff' must be non null.");
  if( *tmp_buff == NULL ) error("poly_improve: invalid tmp buffer.");
 
  *size_local_buff = 0;
  log_nfa = poly_nfa( angles, used, poly, *local_buff, size_local_buff,
                      logNT_poly );
     
  /* try to improve the polygon, by reducing width; first reduce wmax */ 
  for( n=0; n<10; n++ )
    {      
      if( poly->wmax - delta > 0.0 )
        {
          poly->wmax -= delta;
          size_tmp_buff = 0;
          log_nfa_new = poly_nfa( angles, used, poly, *tmp_buff, &size_tmp_buff,
                                  logNT_poly);
           if( log_nfa_new > log_nfa )
            {
              log_nfa = log_nfa_new;
              /* switch addresses of the pointers to the local and tmp buffer */
              swap_ptr_pts( local_buff, tmp_buff );
              (*size_local_buff) = size_tmp_buff;
            }    
        }
    }

  /* now try to reduce dmin */
  for( n=0; n<10; n++ )
    {      
      if( poly->wmin + delta < 0.0 )
        {
          poly->wmin += delta;
          size_tmp_buff = 0;
          log_nfa_new = poly_nfa( angles, used, poly, *tmp_buff, &size_tmp_buff,
                                  logNT_poly);
           if( log_nfa_new > log_nfa )
            {
              log_nfa = log_nfa_new;
              /* switch addresses of the pointers to the local and tmp buffer */
              swap_ptr_pts( local_buff, tmp_buff );
              (*size_local_buff) = size_tmp_buff;
            }    
        }
    }
  return log_nfa;
}


/*----------------------------------------------------------------------------*/
/** Given a polygon as a list of rectangles, compute best validation score.
 */
static double valid_poly( PImageDouble angles, PImageInt used, PolyRect *poly,
                          Point **local_buff, int *size_local_buff, 
                          Point **tmp_buff, double logNT_poly )
{
  double nfa_poly;

  /* check parameters */
  if( (angles == NULL) || (angles->data == NULL) )
    error("valid_poly: invalid input 'angles' image.");
  if( (used == NULL) || (used->data == NULL) )
    error("valid_poly: invalid input 'used' image.");
  if( rec == NULL ) error("valid_poly: invalid input rectangle.");
  if( *local_buff == NULL ) error("valid_poly: invalid local buffer.");
  if( size_local_buff == NULL ) 
    error("valid_poly: 'size_local_buff' must be non null.");
  if( *tmp_buff == NULL ) error("valid_poly: invalid tmp buffer.");
  
  /* compute best validation score */
  nfa_poly = poly_improve( angles, used, poly, local_buff, size_local_buff,
                           tmp_buff, logNT_poly );
   
  /* translate by 0.5 the ends of the segments, so that they correspond to
     pixels' centers */
  for( i=0; i<poly->dim; i++ )
    { 
      poly->rectlist[i].x1 += 0.5; poly->rectlist[i].y1 += 0.5;
      poly->rectlist[i].x2 += 0.5; poly->rectlist[i].y2 += 0.5;
    }  

  return nfa_poly;
}


/*----------------------------------------------------------------------------*/
/*--------- FUNCTIONS FOR CIRCLE RING ESTIMATION AND VALIDATION --------------*/
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/** Compute the widths of a circle ring covering a set of pixels, given the
    centre of the circle and the radius.
 */
static void circ_ring_width( Point *reg, int reg_size, Ring *cring )
{
  int i;
  double w;
  double wmin = -DBL_MAX;
  double wmax = DBL_MAX;

  /* check parameters */
  if( reg == NULL ) error("circ_ring_width: invalid region.");
  if( cring == NULL ) error("circ_ring_width: invalid ring.");
  
  /* compute widths */ 
  for( i=0; i<reg_size; i++ )
    {  
      w = dist( reg[i].x, reg[i].y, cring->cx, cring->cy) - cring->ax;    
      if( w<wmin ) wmin = w;
      if( w>wmax ) wmax = w;
    }
  cring->wmin = wmin; cring->wmax = wmax; 
}


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
  double d;
  double theta0, theta;

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
  *size_buff = 1;
  buff[0].x = x;
  buff[0].y = y;
  used->data[ buff[0].y * used->xsize + buff[0].x ] = label;

  for( i=0; i < *size_buff; i++ )
    for( xx=buff[i].x-1; xx<=buff[i].x+1; xx++ )
      for( yy=buff[i].y-1; yy<=buff[i].y+1; yy++ )
        if( in_image( xx, yy, used ) )
          {
            adr = xx + yy * used->xsize;
            /* check if the point has already been used or visited */
            if( (used->data[adr] != USED) && (used->data[adr] != label) ) 
              {
                /* check if point inside the circle ring by computing the 
                   distance from point to circle and angle */
                d = dist( xx, yy, cring->cx, cring->cy) - cring->ax;
                theta0 = atan2( (double)yy - cring->cy, (double)xx - cring->cx);
                /* make sure angle in [0, 2pi] */
                if( theta0<0 ) theta = theta0 + M_2__PI;
                else theta = theta0;
                
                if( is_in_ring( cring, d, theta ) )
                  {
                    /* check if gradient converges or diverges when computing 
                       the ideal angle of a point on the circle */
                    if( grad_dir == 0 )
                      {
		        if( theta0 > 0 ) theta = -(M_PI-theta0);
		        else theta = M_PI + theta0;			
	              }
                    else theta = theta0;

                    (*err) += norm_angle_diff( angles->data[adr], theta );
                    used->data[adr] = label;
                    buff[*size_buff].x = xx;
                    buff[*size_buff].y = yy;
                    ++(*size_buff); 
                  }
              }
          }          
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
  int n;
  int size_tmp_buff;
  Point *tmp;

  /* check parameters */
  if( (angles == NULL) || (angles->data == NULL) )
    error("circ_ring_improve: invalid input 'angles' image.");
  if( (used == NULL) || (used->data == NULL) )
    error("circ_ring_improve: invalid input 'used' image.");
  if( cring == NULL ) error("circ_ring_improve: invalid input ring.");
  if( buff == NULL ) error("circ_ring_improve: invalid local buffer.");
  if( size_buff == NULL ) 
    error("circ_ring_improve: 'size_local_buff' must be non null.");
  if( *label <= 1 ) error("circ_ring_improve: forbidden label value.");

  (*label)++;
  *size_local_buff = 0;
  log_nfa = circ_ring_nfa( angles, used, cring, x, y, grad_dir, *local_buff, 
                           size_local_buff, *label, logNT_ell );
     
  /* try to improve the ring, by reducing width; first reduce wmax */
  Ring r; 
  copy_ring( cring, &r );
  for( n=0; n<10; n++ )
    {        
      if( r.wmax-delta > 0.0 )
        {
          r.wmax -= delta;
          size_tmp_buff = 0;
          (*label)++;
          log_nfa_new = circ_ring_nfa( angles, used, cring, x, y, grad_dir,
                                       *tmp_buff, &size_local_buff, *label,
                                       logNT_ell );
          if( log_nfa_new > log_nfa )
            {
              cring->wmax -= delta;
              log_nfa = log_nfa_new;
              /* switch addresses of the pointers to the local and tmp buffer */
              swap_ptr_pts( local_buff, tmp_buff );
              (*size_local_buff) = size_tmp_buff;
            }            
        }
    }

  /* now try to reduce wmin */
  ring_copy( cring, &r );
  for( n=0; n<10; n++ )
    {         
      if( r.wmin + delta < 0.0 )
        {
          r.wmin += delta;
          size_tmp_buff = 0;
          (*label)++;
          log_nfa_new = circ_ring_nfa( angles, used, cring, x, y, grad_dir,
                                       *tmp_buff, &size_local_buff, *label,
                                       logNT_ell );
          if( log_nfa_new > log_nfa )
            {
              cring->wmin += delta;
              log_nfa = log_nfa_new;
              /* switch addresses of the pointers to the local and tmp buffer */
              swap_ptr_pts( local_buff, tmp_buff );
              (*size_local_buff) = size_tmp_buff;
            }            
        }  
    }  
  return log_nfa;
}


/*----------------------------------------------------------------------------*/
/** Given a list of points and a fitted circle, compute the circle ring
    covering the points and find best validation score.
 */
static double valid_circ( PImageDouble angles, PImageInt used, Ring *cring,
                          Point *reg, int reg_size, Point **local_buff, 
                          int *size_local_buff, Point **tmp_buff, int dir, 
                          int *label, double logNT_ell )
{
  double nfa_circ;
 
   /* check parameters */
  if( (angles == NULL) || (angles->data == NULL) )
    error("valid_circ: invalid input 'angles' image.");
  if( (used == NULL) || (used->data == NULL) )
    error("valid_circ: invalid input 'used' image.");
  if( cring == NULL ) error("valid_circ: invalid input ring.");
  if( reg == NULL ) error("valid_circ: invalid input region.");
  if( reg_size <= 0 ) error("valid_circ: invalid region size.");
  if( *local_buff == NULL ) error("valid_circ: invalid local buffer.");
  if( size_local_buff == NULL ) 
    error("valid_circ: 'size_local_buff' must be non null.");
  if( *tmp_buff == NULL ) error("valid_circ: invalid tmp buffer.");
  if( *label <= 1 ) error("circ_ring_improve: forbidden label value.");

  /* get the circular ring covering the region of pixels */
  circ_ring_width( reg, reg_size, cring );
 
  /* compute best validation score */
  nfa_circ = circ_ring_improve( angles, used, cring, reg[0].x, reg[0].y, 
                                grad_dir, local_buff, size_local_buff, 
                                tmp_buff, label, logNT_ell);
 
  return nfa_circ;
}


/*----------------------------------------------------------------------------*/
/*--------- FUNCTIONS FOR ELLIPSE RING ESTIMATION AND VALIDATION -------------*/
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/** Compute the widths of an ellipse ring covering a set of pixels, given the
    parameters of the ellipse.
 */
static void ell_ring_width( Point *reg, int reg_size, Ring *ering )
{
  int i;
  double w;
  double wmin = -DBL_MAX;
  double wmax = DBL_MAX;

  /* check parameters */
  if( reg == NULL ) error("ell_ring_width: invalid region.");
  if( ering == NULL ) error("ell_ring_width: invalid ring.");
  
  /* compute widths */ 
  for( i=0; i<reg_size; i++ )
    {  
      w = d_rosin( ering, (double)reg[i].x, (double)reg[i].y );
      if( w<wmin ) wmin = w;
      if( w>wmax ) wmax = w;
    }
  /////////////////MOVE ELSEWHERE if (dmax-dmin > ering->bx) return 0;
  ering->wmin = wmin; ering->wmax = wmax; 
}


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
  int xx, yy, i;
  double d;
  double theta0, theta;

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
  *size_buff = 1;
  buff[0].x = x;
  buff[0].y = y;
  used->data[ buff[0].y * used->xsize + buff[0].x ] = label;

  for( i=0; i < *size_buff; i++ )
    for( xx=buff[i].x-1; xx<=buff[i].x+1; xx++ )
      for( yy=buff[i].y-1; yy<=buff[i].y+1; yy++ )
        if( in_image( xx, yy, used ) )
          {
            adr = xx + yy * used->xsize;
            /* check if the point has already been used or visited */
            if( (used->data[adr] != USED) && (used->data[adr] != label) ) 
              {
                /* check if point inside the ellipse ring by computing the 
                   distance from point to ellipse and angle */
                d = d_rosin( ering, (double)xx, (double)yy );
                theta0 = atan2( (double)yy - ering->cy, (double)xx - ering->cx);
                /* make sure angle in [0, 2pi] */
                if( theta0 < 0 ) theta = theta0 + M_2__PI;
                else theta = theta0;
                
                if( is_in_ring( ering, d, theta ) )
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
                   label, pts, err );

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
  int n;
  int size_tmp_buff;
  Point *tmp;

  /* check parameters */
  if( (angles == NULL) || (angles->data == NULL) )
    error("ell_ring_improve: invalid input 'angles' image.");
  if( (used == NULL) || (used->data == NULL) )
    error("ell_ring_improve: invalid input 'used' image.");
  if( ering == NULL ) error("ell_ring_improve: invalid input ring.");
  if( foci == NULL ) error("ell_ring_improve: invalid input foci.");
  if( buff == NULL ) error("ell_ring_improve: invalid local buffer.");
  if( size_buff == NULL ) 
    error("ell_ring_improve: 'size_local_buff' must be non null.");
  if( *label <= 1 ) error("ell_ring_improve: forbidden label value.");

  (*label)++;
  *size_local_buff = 0;
  log_nfa = ell_ring_nfa( angles, used, ering, foci, x, y, grad_dir, 
                          *local_buff, size_local_buff, *label, logNT_ell );
     
  /* try to improve the ring, by reducing width; first reduce wmax */
  Ring r; 
  copy_ring( ering, &r );
  for( n=0; n<10; n++ )
    {        
      if( r.wmax-delta > 0.0 )
        {
          r.wmax -= delta;
          size_tmp_buff = 0;
          (*label)++;
          log_nfa_new = ell_ring_nfa( angles, used, ering, foci, x, y, grad_dir,
                                      *tmp_buff, &size_local_buff, *label,
                                      logNT_ell );
          if( log_nfa_new > log_nfa )
            {
              ering->wmax -= delta;
              log_nfa = log_nfa_new;
              /* switch addresses of the pointers to the local and tmp buffer */
              swap_ptr_pts( local_buff, tmp_buff );
              (*size_local_buff) = size_tmp_buff;
            }            
        }
    }

  /* now try to reduce wmin */
  ring_copy( ering, &r );
  for( n=0; n<10; n++ )
    {         
      if( r.wmin + delta < 0.0 )
        {
          r.wmin += delta;
          size_tmp_buff = 0;
          (*label)++;
          log_nfa_new = ell_ring_nfa( angles, used, ering, foci, x, y, grad_dir,
                                      *tmp_buff, &size_local_buff, *label,
                                      logNT_ell );
          if( log_nfa_new > log_nfa )
            {
              ering->wmin += delta;
              log_nfa = log_nfa_new;
              /* switch addresses of the pointers to the local and tmp buffer */
              swap_ptr_pts( local_buff, tmp_buff );
              (*size_local_buff) = size_tmp_buff;
            }            
        }  
    }  
  return log_nfa;
}   
 

/*----------------------------------------------------------------------------*/
/** Given a list of points and a fitted ellipse, compute the ellipse ring
    covering the points and find best validation score.
 */
static double valid_ell( PImageDouble angles, PImageInt used, Ring *ering,
                         Point *reg, int reg_size, Point **local_buff, 
                         int *size_local_buff, Point **tmp_buff, int dir, 
                         int *label, double logNT_ell )
{
  double nfa_ell;
  double foci[4]; /* xf1, yf1, xf2, yf2 */
 
  /* check parameters */
  if( (angles == NULL) || (angles->data == NULL) )
    error("valid_ell: invalid input 'angles' image.");
  if( (used == NULL) || (used->data == NULL) )
    error("valid_ell: invalid input 'used' image.");
  if( ering == NULL ) error("valid_circ: invalid input ring.");
  if( reg == NULL ) error("valid_circ: invalid input region.");
  if( reg_size <= 0 ) error("valid_circ: invalid region size.");
  if( *local_buff == NULL ) error("valid_circ: invalid local buffer.");
  if( size_local_buff == NULL ) 
    error("valid_circ: 'size_local_buff' must be non null.");
  if( *tmp_buff == NULL ) error("valid_circ: invalid tmp buffer.");
  if( *label <= 1 ) error("circ_ring_improve: forbidden label value.");

  /* compute foci of the given ellipse */
  ellipse_foci( ering, foci );

  /* get the elliptical ring covering the region of pixels */
  ell_ring_width( reg, reg_size, ering );
 
  /* compute best validation score */
  nfa_ell = ell_ring_improve( angles, used, ering, foci, reg[0].x, reg[0].y, 
                              grad_dir, local_buff, size_local_buff, tmp_buff, 
                              label, logNT_ell);
  return nfa_ell;
}


/*----------------------------------------------------------------------------*/
/*--------- FUNCTIONS FOR GATHERING CONNECTED CONVEX REGIONS OF PIXELS--------*/
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/** Compute the seed for new region growing. All the pixels situated at the end
    of the previously detected rectangle and oriented similarly with it, 
    are considered as seed buffer for new region growing.  
 */
static void px_seed( PImageDouble angles, PImageInt used, Rectangle *rect,  
                     Point *reg, int start, int end, int *reg_buff_size,
                     double *pext, double prec, int label )
{
  int xx, yy, i;
  int adr;
  double xtmp = pext[0], ytmp = pext[1];
  double width = rec->width;
  double l1, l2, l;

  /* check parameters */
  if( (angles == NULL) || (angles->data == NULL) )
    error("px_seed: invalid input 'angles' image.");
  if( (used == NULL) || (used->data == NULL) )
    error("px_seed: invalid input 'used' image.");
  if( rect == NULL ) error("px_seed: invalid input rectangle.");
  if( reg == NULL ) error("px_seed: invalid input region.");
  if( end <= start ) error("px_seed: invalid region size.");

  /* new extreme point will be the end of the rectangle farthest to the previous
     extreme point */
  l1 = fabs( (rec->x1 - xtmp) * rec->dx + (rec->y1 - ytmp) * rec->dy );
  l2 = fabs( (rec->x2 - xtmp) * rec->dx + (rec->y2 - ytmp) * rec->dy );
  if( l1>l2 ) 
    { pext[0] = rec->x1; pext[1] = rec->y1; }
  else
    { pext[0] = rec->x2; pext[1] = rec->y2; }

  /* init size of seed points buffer */
  /* *reg_size_buff = end; */

  /* distance between rectangle centre and extreme point; a new point can be 
     a seed only if it is at a distance farther than this */ 
  l1 = fabs( (pext[0] - rec->cx) * rec->dx + (pext[1] - rec->cy) * rec->dy);

  /* scan previous region to find valid points that are close to the end of 
     the rectangle; the new seed pixels will be among the neighbours of these 
     points. */
  for( i=start; i<end; i++ )
    {
      l = fabs( ( (double)reg[i].x - pext[0] ) * rec.dx + 
                ( (double)reg[i].y - pext[1] ) * rec.dy);

      if( l<width )
        for( xx=reg[i].x-1; xx<=reg[i].x+1; xx++ )
          for( yy=reg[i].y-1; yy<=reg[i].y+1; yy++ )
            {
              /* if point in image */
              if( in_image( xx, yy, used ) )
                {                   
                  adr = yy * used->xsize + xx;                  
                  /* if point not already used or visited or notdef */
                  if( (used->data[adr] != USED) && (used->data[adr] != label ) &&
                      (angles->data[adr] != NOTDEF)
                    {
                      l2 = fabs( (xx - rec->cx) * rec.dx + 
                                 (yy - rec->cy) * rec.dy ); 
                      /* if point's orientation is similar and the point is at
                         an acceptable distance  */
                      if( (angle_diff( angles->data[adr], 
                           rec->theta-M_1_2_PI ) < 3 * prec ) && (l2 >= l1-1) )
                        {
                          reg[ *reg_size_buff ].x = xx;
                          reg[ *reg_size_buff ].y = yy;
                          (*reg_size_buff)++;
                          used->data[adr] = label;
                        }
                    }
                }
            }
    }  

  /* if none of the visited pixels satisfy the conditions, take as seed 
     the end of the rectangle, if it is inside the image */
  if( end == *reg_size_buff )
    {
      xx = round(pext[0]); yy = round(pext[1]);
      if( in_image( xx, yy, used ) )
        {
          reg[ *reg_size_buff ].x = xx; 
          reg[ *reg_size_buff ].y = yy;
          (*reg_size_buff++);
        }
    } 
}


/*----------------------------------------------------------------------------*/
/** Starting from an initial rectangle, scan in one direction to compute a 
    polygonal approximation of a relatively smooth and convex curve.
    A new rectangle is accepted if the angle difference is less than PI/2, 
    and the length ratio is less than 10. These are empirical values.    
 */
static void subcurve( PImageDouble angles, PImageDouble gradmag, PImageInt used, 
                      PolyRect *poly, Point *reg, int start, int *end, 
                      int reg_buff_size, double *pext,  double *spir,
                      double density_th, double prec,int *sgn, int *label, 
                      int *poly_dim_max )
{
  int convex = 1;
  double reg_angle, difang, ang_prev;
  double lenrec, lenth;
  Rectangle rec;
  int i;

  /* check parameters */
  if( (angles == NULL) || (angles->data == NULL) )
    error("subcurve: invalid input 'angles' image.");
  if( (gradmag == NULL) || (gradmag->data == NULL) )
    error("subcurve: invalid input 'gradmag' image.");
  if( (used == NULL) || (used->data == NULL) )
    error("subcurve: invalid input 'used' image.");
  if( rec == NULL ) error("subcurve: invalid input rectangle.");
  if( reg == NULL ) error("subcurve: invalid input region.");
  if( *end <= start ) error("subcurve: invalid region size.");
  if( *label <= 1 ) error("circ_ring_improve: forbidden label value.");


  /* get length and angle of initial rectangle (stored in polygon on first 
     position) */
  lenth = dist( poly->rectlist[0].x1, poly->rectlist[0].y1,
                poly->rectlist[0].x2, poly->rectlist[0].y2 );
  ang_prev = poly->rectlist[0].theta;
 
  /* add new rectangles while the conditions of smoothness and convexity 
     are satisfied */
  while(convex)
    {
      (*label)++;
      start = *end;
      *end = reg_buff_size;
      /* start region grow at the end of the previous rectangle */
      region_grow( angles, used, reg, start, reg_buff_size, end, *label, prec,
                   &reg_angle );       
      if( *end - start > 1 )
        {
          /* estimate and refine rectangle */
          region2rect( gradmag, reg, start, *end, reg_angle, &rec, prec );
          (*label)++;
          if( !refine( gradmag, angles, used, reg, start, end, reg_buff_size,
                       *label, &rec, density_th, prec) ) 
            {
              /* if operation not successful, clean visited points and return */
              mark_img_pts( used, reg, start, *end, NOTUSED );
              *end = start;
              return;
            }

          /* if consecutive rectangles don't have comparable lengths  
             (condition for smoothness): clean and return */
          lenrec = dist(rec.x1,rec.y1,rec.x2,rec.y2);
          if( lenrec < lenth/5.0 ) 
            {
              mark_img_pts( used, reg, start, *end, NOTUSED );
              *end = start;
              return;
            }
          lenth = lenrec;

          /* Convexity rule imposes that consecutive rectangles turn in the same
             direction, hence must have the same sign for angle difference 
             between consecutive rectangles' orientation. If set, verify sign
             before adding a new rectangle. If not set (at first call), 
             then set sign. */
          if( *sgn == 0 ) 
            *sgn = sign( angle_diff_signed( rec.theta, ang_prev ) ); 

          /* get angle difference between current rectangle and previous 
             rectangle */  
          difang = angle_diff_signed( rec.theta, ang_prev );

          /* convexity & non-spiral check */
          if( (sign(difang) == *sgn) && (fabs(difang) < M_1_2_PI) &&
              ( (*spir = *spir + fabs(difang)) <= M_2__PI) )
            { 
              /* if all conditions are met, find next seed points, 
                 and add rectangle to polygon */
              (*label)++;
              reg_buff_size = *end;
              px_seed( angles, used, &rec, reg, start, *end, &reg_buff_size,
                       pext, prec, *label );

              /* add rectangle to polygon */
              poly->dim++;
              add_rect_to_poly( poly, &rec, *poly_dim_max );
              ang_prev = rec.theta; 
            }      
          else /* not convex */
            {
              convex = 0;
              /* substract the contribution of the last rectangle and clean */
              mark_img_pts( used, reg, start, *end, NOTUSED );
              *end = start;
            }
        } /* !if( *end - start > 1 ) */
      else 
        { 
          convex = 0;
          *end = start;
        }     
    }  /* !while(convex) */
}


 


/*----------------------------------------------------------------------------*/
/** Compute a polygonal approximation of a relatively smooth and convex curve, 
    scanning in both directions starting from an initial rectangle.    
 */
static int curve_grow( PDoubleImage gradmag, PImageDouble angles, 
                       PImageInt used, Point *reg, int *reg_size, double **buff, 
                       int *size_buff_max, double density_th, double prec, 
                       int *label, int min_size_ell ) 

{
  double pext[2];
  int sgn = 0;
  double spir = 0.0;
  Rectangle rec;
  double vgg[9];
  int size_first_reg = 0;
  int start, idx_buff;
  double reg_angle;
  int poly_dim_max;

  /* check parameters */
  if( (gradmag == NULL) || (gradmag->data == NULL) )
    error("curve_grow: invalid input 'gradmag' image.");
  if( (angles == NULL) || (angles->data == NULL) )
    error("curve_grow: invalid input 'angles' image.");
  if( (used == NULL) || (used->data == NULL) )
    error("curve_grow: invalid input 'used' image.");
  if( reg == NULL ) error("curve_grow: invalid input region.");
  if( *reg_size <=0 ) error("curve_grow: list must contain at least 1 point.");
  if( *buff == NULL ) error("curve_grow: invalid buffer.");
  if( *label <= 1 ) error("circ_ring_improve: forbidden label value.");
  
  /* init indexes */
  start = 0;
  idx_buff = *reg_size;
  (*label)++;

  /* perform initial region growing, to estimate first rectangle */
  region_grow( angles, used, reg, start, idx_buff, *reg_size, *label, prec, 
               &reg_angle ); 
  
  if( *reg_size <= 2 ) return 0;
  
  /* estimate initial rectangle */
  region2rect( gradmag, reg, start, *reg_size, reg_angle, &rec, prec );

  if( !refine( gradmag, angles, used, reg, start, reg_size, idx_buff, &label,
               &rec, density_th, prec )
    {
      mark_img_pts( used, reg, start, *reg_size, NOTUSED );
      return 0;
    }
  
  /* save size of first region; useful to compute second seed, after scanning
     first direction */
  size_first_reg = *reg_size;
  /* add rectangle to poly */
  poly->dim++; 
  /* initial max number of rectangles in polygon is zero */
  poly_dim_max = 0;
  add_rect_to_poly( poly, &rec, &poly_dim_max );

  /* compute first seed with the first end of the initial rectangle */
  pext[0] = rec.x1; 
  pext[1] = rec.y1;
  (*label)++;
  px_seed( angles, used, &rec, reg, start, *reg_size, &idx_buff, pext, prec,
           *label ); 

  /* lauch scan in first direction */ 
  subcurve( angles, gradmag, used, poly, reg, start, *reg_size, idx_buff, 
            pext, &spir, density_th, prec, &sgn, label, poly_dim_max ); 
  
  /* scan the other direction only if a complete tour was not scanned; this
     condition is to avoid spiral curves */ 
  if( spir < M_2__PI )
    {
      /* compute second seed with the second end of the initial rectangle */  
      pext[0] = rec.x2; 
      pext[1] = rec.y2; 
      idx_buff = *reg_size;
      (*label)++;
      px_seed( angles, used, &rec, reg, start, size_first_reg, &idx_buff, 
               pext, prec, *label );
      /* when scanning the other direction, the sign must be inverted */ 
      sgn = -sgn;
      /* lauch scan in second direction */ 
      subcurve( angles, gradmag, used, poly, reg, start, *reg_size, idx_buff, 
            pext, &spir, density_th, prec, &sgn, label, poly_dim_max );
    }  
  
  return 1;
}




  if (reg_size > min_size_ell) 
    {
      fit_equations(vgg,reg,reg_size);
      /* perform first ellipse fitting, and then circle fitting because 
         circle fitting modifies the equations 'eq' */
      fitellipse(vgg,ering,reg_size);
      fitcircle(vgg,cring,reg_size); 
      *dirc = int_ext_c(cring->cx,cring->cy,reg[0].x,reg[0].y);
      double foci[4];
      check_ellipse(ering);
      ellipse_foci(ering,foci);
      *dire = int_ext_e(reg[0].x,reg[0].y,foci);
      if (*dirc)
        {
          cring->x1 = pext1[0]; cring->y1 = pext1[1];
          cring->x2 = pext2[0]; cring->y2 = pext2[1];
        }
      else
        {
          cring->x1 = pext2[0]; cring->y1 = pext2[1];
          cring->x2 = pext1[0]; cring->y2 = pext1[1];
        }

      if (*dire)
        {
          ering->x1 = pext1[0]; ering->y1 = pext1[1];
          ering->x2 = pext2[0]; ering->y2 = pext2[1];
        }
      else
        {
          ering->x1 = pext2[0]; ering->y1 = pext2[1];
          ering->x2 = pext1[0]; ering->y2 = pext1[1];
        }
      cring->full = 0;
      ering->full = 0;
      cring->ang_start = atan2((double)cring->y1 - cring->cy,
                               (double)cring->x1 - cring->cx);
      if (cring->ang_start < 0) cring->ang_start += M_2__PI; 

      cring->ang_end = atan2((double)cring->y2 - cring->cy,
                             (double)cring->x2 - cring->cx);
      if (cring->ang_end < 0) cring->ang_end += M_2__PI;

      ering->ang_start = atan2((double)ering->y1 - ering->cy,
                               (double)ering->x1 - ering->cx);
      if (ering->ang_start < 0) ering->ang_start += M_2__PI; 

      ering->ang_end = atan2((double)ering->y2 - ering->cy,
                             (double)ering->x2 - ering->cx);
      if (ering->ang_end < 0) ering->ang_end += M_2__PI;
      if (cring->ang_start > cring->ang_end) 
        tang = cring->ang_end + M_2__PI;
      else tang = cring->ang_end;
      
      if (spir > M_PI && tang - cring->ang_start < M_1_2_PI/2.0) 
        swap_ring(cring);
      if (ering->ang_start > ering->ang_end) 
        tang = ering->ang_end + M_2__PI;
      else tang = ering->ang_end;
      if (spir > M_PI && tang - ering->ang_start < M_1_2_PI/2.0) 
        swap_ring(ering);
    }
  return 1;
}


/*----------------------------------------------------------------------------*/
/** Entry point in detector's code. Gets as input an image, and returns two 
    lists of geometric shapes: polygons and ellipses (circles are considered 
    particular ellipses, and segments are particular polygons). 
 */
void ELSDc( PImageDouble in, int *ell_count, Ring *ell, int *poly_count, 
            PImageInt out )
{
  double ang_th = 22.5;     /* gradient angle tolerance in degrees */
  double prec;              /* radian precision */
  double p;
  double rho = 5.0;         /* pixels with gradient magnitude under this 
                               threshold are not considered, as their (weak) 
                               orientation might be too affected by noise  */
  int grad_dir;             /* store direction of the gradient wrt center of 
                               circle or ellipse */
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
  PImageDouble imgauss;     /* smooth version of the original image, used during
                               candidate selection */
  Rectangle rec;            /* rectangle (line segment) parameters  */
  PolyRect poly;            /* parameters of polygon as list of rectangles */
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
  int size_new_buff;        /* buffer to store the points of a new primitive; 
                               if the new primitive is more meaningful than 
                               the most meaningful primitive found so far, 
                               deem new primitive as most meaningful by 
                               switching the addresses with best_buff   */
  Point *tmp_buff;          /* buffer to store the points of a new version of a
                               primitive, computed during improve operation; 
                               if the new version is more meaningful than 
                               the most meaningful version found so far, 
                               deem new version as most meaningful by 
                               switching the addresses with new_buff   */


  /* allocate and init variables */
  prec = M_PI*ang_th/180.0;
  p = ang_th/180.0;
  xsize = angles->xsize;
  ysize = angles->ysize;
  poly->dim = 0;  
  poly->rectlist = NULL;
  
  buff_fit = (double*) malloc ( sizeof (double) );
  size_buff_fit = 1;

  reg = (Point *) malloc( xsize * ysize * sizeof(Point) );
  if( reg == NULL ) error("ELSDc: not enough memory.");
  best_buff = (Point *) malloc( xsize * ysize * sizeof(Point) );
  if( best_buff == NULL ) error("ELSDc: not enough memory.");
  new_buff = (Point *) malloc( xsize * ysize * sizeof(Point) );
  if( new_buff == NULL ) error("ELSDc: not enough memory.");
  tmp_buff = (Point *) malloc( xsize * ysize * sizeof(Point) );
  if( tmp_buff == NULL ) error("ELSDc: not enough memory.");

  used = new_image_int_ini(xsize,ysize,NOTUSED);

  /* Number of tests for ellipse arcs. We must count (3*xsize)*(3*ysize) 
     possible values for the center, as center can be outside the image too,
     (xsize*ysize) possible values for the two axes, sqrt(xsize*ysize) possible
     values for orientation, sqrt(xsize*ysize) possible values for width, 
     (xsize*ysize) possible values for delimiting angles. This makes 
     9*(xsize*ysize)^4. Since we have polygons and ellipses, i.e. two families,
     we must multiply this by 2. So NT_ell = 2*9*(xsize*ysize)^4. 
     In log, we get:   */
  logNT_ell = 4.0 * ( log10( (double)xsize ) + log10( (double)ysize ) ) + 
                      log10(9.0) + log10(2.0); 

  logNT_seg = 5.0 * ( log10((double)xsize)+log10((double)ysize))/2.0 + 
                    log10(2.0); /* N^5 */

  min_size_ell =round( (-logNT_ell-mlog10eps) / log10(p) );
  min_size_seg =(int)((-logNT_seg+log10(eps))/log10(p));

  /* perform gaussian smoothing */
  imgauss = gaussian_sampler(in,1.0,0.6);
  
  /* compute gradient magnitude and orientation  for the smooth image */
  ll_angle_full(imgauss,rho,&list_p,&mem_p,n_bins,max_grad);

  /* compute gradient magnitude and orientation  for the original image */
  ll_angle(image,rho);

  
  msize = 100*max(xsize,ysize);

  /* begin primitive detection */
  for(;list_p; list_p = list_p->next)
    { 
      if(used->data[list_p->y*used->xsize+list_p->x] != USED &&
         angles->data[list_p->y*angles->xsize+list_p->x] != NOTDEF)
	{
          /* init some variables */ 	
          nfa_poly = nfa_circ = nfa_ell = nfa_seg = mlog10eps; 
          reg_size = reg_size_buff = 1; reg_size0 = 0;
          USEDTMP++;
          size_poly = size_poly0 = size_circ = size_ell = size_seg = 0; 
          poly->dim = 0; 
          reg[0].x = list_p->x; reg[0].y = list_p->y;
          /* BUILD HYPOTHESES: polygon, circle, ellipse */
          if (!curve_grow(prec,p,min_size_ell,min_size_seg,logNT_seg,
              mlog10eps,density_th,poly,cring,ering,&dirc,&dire))
            continue;
          poly_profile(poly);

          /* VALIDATE HYPOTHESES */
          //printf("-------- Validation -------\n");
          logNT_poly = (poly->dim + 1+1.0/2.0) * (log10((double)used->xsize) + 
               log10((double)used->ysize)) + (poly->dim+1.0)*log10(2.0);
          min_size_poly =(int)((-logNT_poly - mlog10eps)/log10(p));
          if (reg_size>min_size_poly)
            nfa_poly = valid_poly(poly,p,logNT_poly,min_size_poly,mlog10eps);

          size_circ = size_ell = size_seg = size_poly;
          if (reg_size>min_size_ell)
            {
              if (check_circle(cring))
                nfa_circ = valid_circ(cring,dirc,prec,p,logNT_ell,min_size_ell,
                                    mlog10eps,density_th,msize);
              size_ell = size_seg = size_circ;
              if (check_ellipse(ering))
                nfa_ell = valid_ell(ering,dire,prec,p,logNT_ell,min_size_ell,
                                    mlog10eps,density_th,msize);
            }
          size_seg = size_ell;
          if (reg_size > min_size_seg)
            nfa_seg = valid_seg(seg,p,prec,logNT_seg,min_size_seg,mlog10eps);
          
          //printf("-------- NFA values -------\n");
          //printf("Seg %f Poly %f Cir %f Ell %f \n",nfa_seg,nfa_poly,nfa_circ,
                // nfa_ell);
          //printf("%d %d %d %d \n",size_poly,size_circ,size_ell,size_seg);    
          //printf("%d\n",USEDTMP);
          /* MODEL SELECTION BY COMPARING NFAs */
/*ellipse*/if(nfa_ell > mlog10eps && nfa_ell > nfa_poly && 
              nfa_ell > nfa_circ && nfa_ell > nfa_seg) 
             {
	       //printf("ellipse\n");
               (*ell_count)++;//write_svg_poly(svg,poly);
               write_svg_ellipse(fe,svg,ering);
               mark_used(size_circ,size_ell);
             }
/* circle */else if(nfa_circ > mlog10eps && nfa_circ > nfa_poly &&
                    nfa_circ >= nfa_ell && nfa_circ > nfa_seg) 
              { 
                //printf("ellipse\n");
                (*ell_count)++;//write_svg_poly(svg,poly);
                write_svg_circle(fe,svg,cring);
                mark_used(size_poly,size_circ); 
              }
/* polygon */  else if(nfa_poly > mlog10eps && nfa_poly > nfa_circ && 
                    nfa_poly > nfa_ell && nfa_poly > nfa_seg) 
              {
                //printf("polygon\n");
                (*poly_count)++;
		write_svg_poly(svg,poly);             
                mark_used(0,size_poly);
              }
            else if(nfa_seg > mlog10eps && nfa_seg > nfa_circ && 
                    nfa_seg > nfa_ell && nfa_seg >= nfa_poly) 
              {
                //printf("segment\n");
                (*seg_count)++;
                seg->x1 += 0.5; seg->y1 += 0.5;
                seg->x2 += 0.5; seg->y2 += 0.5;
		write_svg_seg(svg,seg);
                //write_svg_poly(svg,poly);                    
                mark_used(size_ell,size_seg);
              }
        }/* IF USED */
    }/* FOR LIST */		    
  free_image_double(imgauss);
  free_image_double(image);
  free_image_double(gradx); free_image_double(grady);
  free_image_double(gradmag); free_image_double(angles);
  free_image_double(angles_orig);  
  free_image_int(used);
  free(reg); free(pts_used);  
  free(poly->rectlist); free(poly); free(ering); free(cring); free(seg); 
  free(gBufferDouble); free(gBufferInt);
  free(mem_p);
  fclose(fe);
  fclose_svg(svg);
}
