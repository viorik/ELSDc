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
#define LOG10_1_2_M_2__PI  0.39908993418
#define LOG10_EXP1         0.43429448190
static double NFAc( int n, double k, double logNT, FILE *fdebug )
{
  fprintf(fdebug,"NFAc %d %lf %lf \n", n, k, logNT);
  /* check parameters */
  if( n < 0) error("NFAc: 'n' must be strictly positive.");
  if( n==0 ) return -logNT;
  double logNFAC = logNT;
  //if( k>0 ) 
  k = max(k,0.01);
  logNFAC = logNT + n * log10(k) - LOG10_1_2_M_2__PI - 
            (n+0.5) * log10(n) + n * LOG10_EXP1 ;

  fprintf(fdebug,"end NFAc %lf \n", -logNFAC);
  return -logNFAC;
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
                         int label_start, double prec, double *reg_angle, 
                         FILE *fdebug )
{
  double sumdx = 0.0, sumdy = 0.0;
  double ang;
  int xx, yy, i;
  int adr;
  fprintf(fdebug,"region grow %d %d %d %d\n",start,idx_buff,*end,label);
  /* check parameters */
  if( (angles == NULL) || (angles->data == NULL) )
    error("region_grow: invalid image 'angles'.");
  if( (used == NULL) || (used->data == NULL) )
    error("region_grow: invalid image 'used'.");
  if( reg == NULL ) error("region_grow: invalid input region.");
  if( start >= *end ) error("region_grow: invalid indexes.");
  if( prec < 0.0 ) error("region_grow: 'prec' must be positive.");
  
  /* Initialise *reg_angle with the value of the pixels 
     orientations in the buffer part of reg. */
  for( i=start; i<idx_buff; i++ )
    {
      xx = reg[i].x; yy = reg[i].y;
      adr = xx + yy*angles->xsize;
      sumdx += cos( angles->data[adr] );
      sumdy += sin( angles->data[adr] );
      used->data[adr] = label;
    }
  *reg_angle = atan2( sumdy, sumdx );
  
  /* try neighbors as new region points */
  for( i=start; i<*end; i++ )
    for( xx=reg[i].x-1; xx<=reg[i].x+1; xx++ )
      for( yy=reg[i].y-1; yy<=reg[i].y+1; yy++ )
        if( in_image( xx, yy, used->xsize, used->ysize ) )
          {
            adr = yy * angles->xsize + xx;
            if( (used->data[adr] != USED) && (used->data[adr] != label) && 
                (used->data[adr] < label_start) )
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
                    fprintf(fdebug,"%d %d\n",xx,yy);
                  }
         
              }
          }
}




/*----------------------------------------------------------------------------*/
/*--------- FUNCTIONS FOR GATHERING CONNECTED CONVEX REGIONS OF PIXELS -------*/
/*----------------------------------------------------------------------------*/


/*------------------ FUNCTIONS FOR RECTANGLE ESTIMATION ----------------------*/


/*----------------------------------------------------------------------------*/
/** Estimate the angle of the rectangle covering the pixels listed in 'reg' 
    between 'start' and 'end' indexes. 'cx' and 'cy' are the centroid 
    coordinates of the point set. The modulus of the gradient magnitude 
    contained in 'gradmag' is used to weight pixels' contributions. 
 */
static double get_theta( PImageDouble gradmag, Point *reg, int start, int end, 
                         double cx, double cy, double reg_angle, double prec, 
                         FILE *fdebug )
{
  double lambda1, lambda2, tmp, theta, weight, sum;
  double Ixx = 0.0;
  double Iyy = 0.0;
  double Ixy = 0.0;
  int i;
  fprintf(fdebug,"get_theta start %d end %d \n",start,end);
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
  fprintf(fdebug,"end get-theta %lf \n",theta);
  return theta;
}


/*----------------------------------------------------------------------------*/
/** Compute a rectangle that covers a connected region of points. A rectangle is
    defined by end points, center, orientation, width. Additional fields 
    (precision of approximation) are useful for validation. 
 */
static void region2rect( PImageDouble gradmag, Point *reg, int start, int end,
                         double reg_angle, Rectangle *rec, double prec, 
                         FILE *fdebug )
{
  double cx, cy, dx, dy, l, w, theta;
  double weight,sum,l_min,l_max,w_min,w_max;
  int i, adr;
  fprintf(fdebug,"region2rect start %d end %d \n",start,end);
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
  theta = get_theta( gradmag, reg, start, end, cx, cy, reg_angle, prec, fdebug );

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
  rec->cx = cx;
  rec->cy = cy;
  rec->theta = theta;
  rec->dx = dx;
  rec->dy = dy;
  rec->prec = prec;

  /* correct if rectangle too thin */
  if( rec->width < 1.0 )
    { 
      rec->width = 1.0; rec->wmin = -0.5; rec->wmax = 0.5;
    }
  fprintf(fdebug,"end reg2rect \n");
}


/*----------------------------------------------------------------------------*/
/** Reduce the region size, by eliminating the points far from the starting 
    point, until that leads to rectangle with the right density of region points
    or to discard the region if too small.
 */
static int reduce_region_radius( PImageDouble gradmag, PImageDouble angles,
                                 PImageInt used, Point *reg, int start, 
                                 int *end, double reg_angle, Rectangle *rec, 
                                 int x, int y, double density_th, double prec, 
                                 FILE *fdebug )
{
  double density, rad1, rad2, rad;
  int i;
  fprintf(fdebug,"reduce reg radius start %d end %d \n",start,*end);
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
      region2rect( gradmag, reg, start, *end, reg_angle, rec, prec, fdebug );

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
                   Point *reg, int start, int idx_buff, int *end, int *label,
                   int label_start, Rectangle *rec, double density_th, 
                   double prec, FILE *fdebug )
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
  fprintf(fdebug,"refine start %d buff %d end %d dens %lf\n",start,idx_buff,*end, density);
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
  region_grow( angles, used, reg, start, idx_buff, end, *label, label_start, 
               tau, &reg_angle, fdebug );

  /* if the region is too small, reject */
  if( *end - start <= 2 ) return FALSE;

  /* re-compute rectangle */
  region2rect( gradmag, reg, start, *end, reg_angle, rec, prec, fdebug );

  /* re-compute region points density */
  density = (double) (*end - start) /
                     ( dist(rec->x1,rec->y1,rec->x2,rec->y2) * rec->width);

  /*------ Second try: reduce region radius ------*/
  if( density < density_th )
    return reduce_region_radius( gradmag, angles, used, reg, start, end, 
                                 reg_angle, rec, xc, yc, density_th, prec, fdebug );
  return TRUE;
}


/*----------------------------------------------------------------------------*/
/*------------------ FUNCTIONS FOR POLYGON AND CURVE ESTIMATION --------------*/
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/** Compute the seed for new region growing. All the pixels situated at the end
    of the previously detected rectangle and oriented similarly with it, 
    are considered as seed buffer for new region growing.  
 */
static void px_seed( PImageDouble angles, PImageInt used, Rectangle *rec,  
                     Point *reg, int start, int *idx_buff, int end,
                     double *pext, double prec, int label, FILE *fdebug )
{
  int xx, yy, i;
  int adr;
  double xtmp = pext[0], ytmp = pext[1];
  double width = rec->width;
  double l1, l2, l;
  fprintf(fdebug,"px_seed start %d buff %d end %d \n",start,*idx_buff,end);
  /* check parameters */
  if( (angles == NULL) || (angles->data == NULL) )
    error("px_seed: invalid input 'angles' image.");
  if( (used == NULL) || (used->data == NULL) )
    error("px_seed: invalid input 'used' image.");
  if( rec == NULL ) error("px_seed: invalid input rectangle.");
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

  /* distance between rectangle centre and extreme point; a new point can be 
     a seed only if it is at a distance farther than this */ 
  l1 = fabs( (pext[0] - rec->cx) * rec->dx + (pext[1] - rec->cy) * rec->dy);

  /* scan previous region to find valid points that are close to the end of 
     the rectangle; the new seed pixels will be among the neighbours of these 
     points. */
  for( i=start; i<end; i++ )
    {
      l = fabs( ( (double)reg[i].x - pext[0] ) * rec->dx + 
                ( (double)reg[i].y - pext[1] ) * rec->dy);

      if( l<width )
        for( xx=reg[i].x-1; xx<=reg[i].x+1; xx++ )
          for( yy=reg[i].y-1; yy<=reg[i].y+1; yy++ )
            {
              /* if point in image */
              if( in_image( xx, yy, used->xsize, used->ysize ) )
                {                   
                  adr = yy * used->xsize + xx;                  
                  /* if point not already used or visited or notdef */
                  if( (used->data[adr] != USED) && (used->data[adr] != label ) &&
                      (angles->data[adr] != NOTDEF) )
                    {
                      l2 = fabs( (xx - rec->cx) * rec->dx + 
                                 (yy - rec->cy) * rec->dy ); 
                      /* if point's orientation is similar and the point is at
                         an acceptable distance  */
                      if( (angle_diff( angles->data[adr], 
                           rec->theta-M_1_2_PI ) < 3 * prec ) && (l2 >= l1-1) )
                        {
                          reg[ *idx_buff ].x = xx;
                          reg[ *idx_buff ].y = yy;
                          ++(*idx_buff);
                          used->data[adr] = label;
                        }
                    }
                }
            }
    }  

  /* if none of the visited pixels satisfy the conditions, take as seed 
     the end of the rectangle, if it is inside the image */
  if( end == *idx_buff )
    {
      xx = round(pext[0]); yy = round(pext[1]);
      if( in_image( xx, yy, used->xsize, used->ysize ) )
        {
          reg[ *idx_buff ].x = xx; 
          reg[ *idx_buff ].y = yy;
          ++(*idx_buff);
        }
    } 
}


static void px_seed2( PImageDouble angles, PImageInt used, Rectangle *rec,  
                     Point *reg, int start, int *idx_buff, int end,
                     double *pext, double prec, int label, FILE *fdebug )
{
  int xx, yy, i;
  int adr;
  double xtmp = pext[0], ytmp = pext[1];
  double width = rec->width;
  double l1, l2, l;
  fprintf(fdebug,"px_seed start %d buff %d end %d \n",start,*idx_buff,end);
  /* check parameters */
  if( (angles == NULL) || (angles->data == NULL) )
    error("px_seed: invalid input 'angles' image.");
  if( (used == NULL) || (used->data == NULL) )
    error("px_seed: invalid input 'used' image.");
  if( rec == NULL ) error("px_seed: invalid input rectangle.");
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

  /* distance between rectangle centre and extreme point; a new point can be 
     a seed only if it is at a distance farther than this */ 
  l1 = fabs( (pext[0] - rec->cx) * rec->dx + (pext[1] - rec->cy) * rec->dy);

  /* scan previous region to find valid points that are close to the end of 
     the rectangle; the new seed pixels will be among the neighbours of these 
     points. */
  for( i=start; i<end; i++ )
    {
      l = fabs( ( (double)reg[i].x - pext[0] ) * rec->dx + 
                ( (double)reg[i].y - pext[1] ) * rec->dy);

      if( l<width )
        for( xx=reg[i].x-1; xx<=reg[i].x+1; xx++ )
          for( yy=reg[i].y-1; yy<=reg[i].y+1; yy++ )
            {
              /* if point in image */
              if( in_image( xx, yy, used->xsize, used->ysize ) )
                {                   
                  adr = yy * used->xsize + xx;                  
                  /* if point not already used or visited or notdef */
                  if( (used->data[adr] != USED) && (used->data[adr] != label ) &&
                      (angles->data[adr] != NOTDEF) )
                    {
                      l2 = fabs( (xx - rec->cx) * rec->dx + 
                                 (yy - rec->cy) * rec->dy ); 
                      /* if point's orientation is similar and the point is at
                         an acceptable distance  */
                      if( (angle_diff( angles->data[adr], 
                           rec->theta-M_1_2_PI ) < 3 * prec ) && (l2 >= l1-1) )
                        {
                          reg[ *idx_buff ].x = xx;
                          reg[ *idx_buff ].y = yy;
                          ++(*idx_buff);
                          used->data[adr] = label;
                        }
                    }
                }
            }
    }  

  /* if none of the visited pixels satisfy the conditions, take as seed 
     the end of the rectangle, if it is inside the image */
  if( end == *idx_buff )
    {
      xx = round(pext[0]); yy = round(pext[1]);
      if( in_image( xx, yy, used->xsize, used->ysize ) )
        {
          reg[ *idx_buff ].x = xx; 
          reg[ *idx_buff ].y = yy;
          ++(*idx_buff);
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
                      PolyRect *poly, Point *reg, int start, int idx_buff,
                      int *end, double *pext,  double *spir, double density_th,
                      double prec,int *sgn, int *label, int label_start, 
                      FILE *fdebug )
{
  int convex = 1;
  double reg_angle, difang, ang_prev;
  double lenrec, lenth;
  int buf_tmp = 0;
  Rectangle rec;
  int i;
  fprintf(fdebug,"subcurve start %d buff %d end %d \n",start,idx_buff,*end);
  /* check parameters */
  if( (angles == NULL) || (angles->data == NULL) )
    error("subcurve: invalid input 'angles' image.");
  if( (gradmag == NULL) || (gradmag->data == NULL) )
    error("subcurve: invalid input 'gradmag' image.");
  if( (used == NULL) || (used->data == NULL) )
    error("subcurve: invalid input 'used' image.");
  if( reg == NULL ) error("subcurve: invalid input region.");
  if( *end <= start ) error("subcurve: invalid region size.");
  if( *label <= 1 ) error("subcurve: forbidden label value.");


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
      *end = idx_buff;
      buf_tmp = *end - start;
      /* start region grow at the end of the previous rectangle */
      region_grow( angles, used, reg, start, idx_buff, end, *label, label_start,
                   prec, &reg_angle, fdebug );       
      if( *end - start > buf_tmp )
        {
          /* estimate and refine rectangle */
          region2rect( gradmag, reg, start, *end, reg_angle, &rec, prec, fdebug );
          (*label)++;
          if( !refine( gradmag, angles, used, reg, start, idx_buff, end,
                       label, label_start, &rec, density_th, prec, fdebug) ) 
            {
              /* if operation not successful, clean visited points and return */
              mark_img_pts( used, reg, start, *end, NOTUSED );
              *end = start;
              fprintf(fdebug,"not refine\n");
              break;
            }
          if( *end - start <= buf_tmp ) {fprintf(fdebug,"no add refine\n"); break;}
          fprintf(fdebug,"refined_rectangle\n");
          write_rectangle( fdebug, &rec );

          /* if consecutive rectangles don't have comparable lengths  
             (condition for smoothness): clean and return */
          lenrec = dist(rec.x1,rec.y1,rec.x2,rec.y2);
          if( lenrec < lenth/5.0 ) 
            {
              mark_img_pts( used, reg, start, *end, NOTUSED );
              *end = start;
              fprintf(fdebug,"too short\n");
              break;
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
              idx_buff = *end;
              px_seed( angles, used, &rec, reg, start, &idx_buff, *end,
                       pext, prec, *label, fdebug );

              /* add rectangle to polygon */
              poly->dim++;
              add_rect_to_polyrect( poly, &rec );
              ang_prev = rec.theta; 
            }      
          else /* not convex */
            {
              convex = 0;
              /* substract the contribution of the last rectangle and clean */
              mark_img_pts( used, reg, start, *end, NOTUSED );
              *end = start;
              fprintf(fdebug,"not convex\n");
            }
        } /* !if( *end - start > 1 ) */
      else 
        { 
          fprintf(fdebug,"not convex no add\n");
          convex = 0;
          *end = start;
        }     
    }  /* !while(convex) */
  fprintf(fdebug,"subcurve end\n");
} 


/*----------------------------------------------------------------------------*/
/** Compute a polygonal approximation of a relatively smooth and convex curve, 
    scanning in both directions starting from an initial rectangle.    
 */
static int curve_grow( PImageDouble gradmag, PImageDouble angles, 
                       PImageInt used, Point *reg, int *reg_size, 
                       double density_th, double prec, PolyRect *poly, 
                       int *label, double *pext1, double *pext2, double *spir, 
                       FILE *fdebug ) 
{
  int sgn = 0;
  Rectangle rec;
  double vgg[9];
  int size_first_reg = 0;
  int label_start;
  int start, idx_buff;
  double reg_angle;
  
  /* check parameters */
  if( (gradmag == NULL) || (gradmag->data == NULL) )
    error("curve_grow: invalid input 'gradmag' image.");
  if( (angles == NULL) || (angles->data == NULL) )
    error("curve_grow: invalid input 'angles' image.");
  if( (used == NULL) || (used->data == NULL) )
    error("curve_grow: invalid input 'used' image.");
  if( reg == NULL ) error("curve_grow: invalid input region.");
  if( *reg_size <=0 ) error("curve_grow: list must contain at least 1 point.");
  if( poly == NULL ) error("curve_grow: poly must be non null.");
  if( *label <= 1 ) error("curve_grow: forbidden label value.");
  
  /* init indexes */
  start = 0;
  idx_buff = *reg_size;
  //(*label)++;
  label_start = (*label);
  /* perform initial region growing, to estimate first rectangle */
  region_grow( angles, used, reg, start, idx_buff, reg_size, *label, label_start,
               prec, &reg_angle, fdebug ); 

  if( *reg_size <= 2 ) return 0;
  
  /* estimate initial rectangle */
  region2rect( gradmag, reg, start, *reg_size, reg_angle, &rec, prec, fdebug );
  
  /* refine rectangle if density of aligned points is less than density 
     threshold */
  if( !refine( gradmag, angles, used, reg, start, idx_buff, reg_size, label, 
               label_start, &rec, density_th, prec, fdebug ) )
    {
      mark_img_pts( used, reg, start, *reg_size, NOTUSED );
      return 0;
    }
  fprintf(fdebug,"refined_rectangle\n");
  write_rectangle( fdebug, &rec );

  /* save size of first region; useful to compute second seed, after scanning
     first direction */
  size_first_reg = *reg_size;

  /* add rectangle to poly */
  poly->dim++; 
  add_rect_to_polyrect( poly, &rec );

  /* compute first seed with the first end of the initial rectangle */
  pext2[0] = rec.x1; 
  pext2[1] = rec.y1;
  (*label)++;
  idx_buff = *reg_size;
  px_seed( angles, used, &rec, reg, start, &idx_buff, *reg_size, pext2, prec,
           *label, fdebug ); 

  /* lauch scan in first direction */ 
  subcurve( angles, gradmag, used, poly, reg, start, idx_buff, reg_size, 
            pext2, spir, density_th, prec, &sgn, label, label_start, fdebug ); 
  
  /* scan the other direction only if a complete tour was not scanned; this
     condition is to avoid spiral curves */
  pext1[0] = rec.x2; 
  pext1[1] = rec.y2; 
  if( *spir < M_2__PI )
    {
      /* compute second seed with the second end of the initial rectangle */  
 
      idx_buff = *reg_size;
      (*label)++;
      px_seed( angles, used, &rec, reg, start, &idx_buff, size_first_reg, 
               pext1, prec, *label, fdebug );
      /* when scanning the other direction, the sign must be inverted */ 
      sgn = -sgn;
      /* lauch scan in second direction */ 
      subcurve( angles, gradmag, used, poly, reg, start, idx_buff, reg_size, 
                pext1, spir, density_th, prec, &sgn, label, label_start, fdebug );
    }    
  return 1;
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
                              Point *buff, int *size_buff, FILE *fdebug )
{
  RectIter *i;
  int adr;
  double theta;
  //fprintf(fdebug,"rect_count_error \n");
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
  fprintf(fdebug,"end rect_count_error %d %f \n", *pts, *err );
}


/*----------------------------------------------------------------------------*/
/** Given a list of points, estimate segment. One segment is a particular case 
    of a polygon.
 */
static void get_seg( PImageDouble angles, PImageDouble gradmag, 
                     PImageInt used, PolyRect *poly, Point *reg, 
                     int reg_size, double prec, FILE *fdebug )
{
  double sumdx = 0.0, sumdy = 0.0;
  double reg_angle;
  int i;
  int adr;
  Rectangle rec;
  fprintf(fdebug,"get_seg \n");
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
  region2rect( gradmag, reg, 0, reg_size, reg_angle, &rec, prec, fdebug ); 
  
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
                              PolyRect *poly, int *pts, double *err, 
                              Point *buff, int *size_buff, FILE *fdebug ) 
{
  int pts_rect;
  double err_rect;
  int i;
  fprintf(fdebug,"poly_count_error \n");
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
                        &err_rect, buff, size_buff, fdebug );
      (*pts) += pts_rect;
      (*err) += err_rect;
    }  
}


/*----------------------------------------------------------------------------*/
/** Iteration to improve NFA: count additive error and compute NFA. 
 */
static double poly_nfa( PImageDouble angles, PImageInt used, PolyRect *poly, 
                        Point *buff, int *size_buff, double logNT_poly, 
                        FILE *fdebug )
{
  int pts;
  double err;
  fprintf(fdebug,"poly_nfa \n");
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
  poly_count_error( angles, used, poly, &pts, &err, buff, size_buff, fdebug );
  /* return computed NFA */
  return NFAc( pts, err, logNT_poly, fdebug );
}


/*----------------------------------------------------------------------------*/
/** Try some polygon variations to improve NFAc value.
 */
static double poly_improve( PImageDouble angles, PImageInt used, PolyRect *poly,
                            Point **local_buff, int *size_local_buff, 
                            Point **tmp_buff, double logNT_poly, FILE *fdebug )
{
  double delta = 0.5;
  double log_nfa,log_nfa_new;
  int n;
  int size_tmp_buff;
  Point *tmp;
  fprintf(fdebug,"poly_improve \n");
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
  log_nfa = poly_nfa( angles, used, poly, *local_buff, size_local_buff,
                      logNT_poly, fdebug );
     
  /* try to improve the polygon, by reducing width; first reduce wmax */ 
  for( n=0; n<10; n++ )
    {      
      if( poly->wmax - delta > delta )
        {
          poly->wmax -= delta;
          size_tmp_buff = 0;
          fprintf(fdebug, "wmin %f wmax %f\n",  poly->wmin,  poly->wmax );
          log_nfa_new = poly_nfa( angles, used, poly, *tmp_buff, &size_tmp_buff,
                                  logNT_poly, fdebug);
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
      if( poly->wmin + delta < -delta )
        {
          poly->wmin += delta;
          size_tmp_buff = 0;
           fprintf(fdebug, "wmin %f wmax %f\n",  poly->wmin,  poly->wmax );
          log_nfa_new = poly_nfa( angles, used, poly, *tmp_buff, &size_tmp_buff,
                                  logNT_poly, fdebug);
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
/*------------------- FUNCTIONS FOR CIRCLE RING VALIDATION -------------------*/
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/** Compute additive normalised per-pixel error inside a circular ring, and 
    stack visited points so they can be marked as 'used' if the circle ring is
    the most meaningful primitive.
 */
static void circ_count_error( PImageDouble angles, PImageInt used, Ring *cring,
                              int x, int y, int grad_dir, Point *buff, 
                              int *size_buff, int label, int *pts, double *err, FILE *fdebug )
{
  int xx, yy, i, adr;
  double d;
  double theta0, theta;
  fprintf(fdebug,"circ_count_error \n");
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
        if( in_image( xx, yy, used->xsize, used->ysize ) )
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
  *pts = *size_buff;         
}


/*----------------------------------------------------------------------------*/
/** Iteration to improve NFA: count additive error and compute NFA. 
 */
static double circ_ring_nfa( PImageDouble angles, PImageInt used, Ring *cring,
                             int x, int y, int grad_dir, Point *buff, 
                             int *size_buff, int label, double logNT_ell, FILE *fdebug )
{
  int pts;
  double err;
  fprintf(fdebug,"circ_ring_nfa \n");
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
                    label, &pts, &err, fdebug );

  /* return computed NFA */
  return NFAc( pts, err, logNT_ell, fdebug );
}


/*----------------------------------------------------------------------------*/
/** Try some circle ring variations to improve NFAc value.
 */
static double circ_ring_improve( PImageDouble angles, PImageInt used, 
                                 Ring *cring, int x, int y, int grad_dir, 
                                 Point **local_buff, int *size_local_buff, 
                                 Point **tmp_buff, int *label, double logNT_ell, FILE *fdebug)
{    
  double log_nfa, log_nfa_new;
  double delta = 0.5;
  int n;
  int size_tmp_buff;
  Point *tmp;
  fprintf(fdebug,"circ_ring_improve \n");
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

  (*label)++;
  *size_local_buff = 0;
  log_nfa = circ_ring_nfa( angles, used, cring, x, y, grad_dir, *local_buff, 
                           size_local_buff, *label, logNT_ell, fdebug );
     
  /* try to improve the ring, by reducing width; first reduce wmax */
  Ring r; 
  copy_ring( cring, &r );
  for( n=0; n<10; n++ )
    {        
      if( r.wmax-delta > delta )
        {
          r.wmax -= delta;
          size_tmp_buff = 0;
          (*label)++;
          fprintf(fdebug, "wmin %f wmax %f\n",  r.wmin, r.wmax );
          log_nfa_new = circ_ring_nfa( angles, used, &r, x, y, grad_dir,
                                       *tmp_buff, &size_tmp_buff, *label,
                                       logNT_ell, fdebug );
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
  copy_ring( cring, &r );
  for( n=0; n<10; n++ )
    {         
      if( r.wmin + delta < -delta )
        {
          r.wmin += delta;
          size_tmp_buff = 0;
          (*label)++;
          fprintf(fdebug, "wmin %f wmax %f\n",  r.wmin, r.wmax );
          log_nfa_new = circ_ring_nfa( angles, used, &r, x, y, grad_dir,
                                       *tmp_buff, &size_tmp_buff, *label,
                                       logNT_ell, fdebug );
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
                             double *err, FILE *fdebug )
{
  int xx, yy, i, adr;
  double d;
  double theta0, theta;
  fprintf(fdebug,"ell_count_error \n");
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
        if( in_image( xx, yy, used->xsize, used->ysize ) )
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
  *pts = *size_buff;     
}


/*----------------------------------------------------------------------------*/
/** Iteration to improve NFA: count additive error and compute NFA. 
 */
static double ell_ring_nfa( PImageDouble angles, PImageInt used, Ring *ering,
                            double *foci, int x, int y, int grad_dir, 
                            Point *buff, int *size_buff, int label, 
                            double logNT_ell, FILE *fdebug )
{
  int pts;
  double err;
  fprintf(fdebug,"ell_ring_nfa \n");
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
                   label, &pts, &err, fdebug );

  /* return computed NFA */
  return NFAc( pts, err, logNT_ell, fdebug );
}


/*----------------------------------------------------------------------------*/
/** Try some ellipse ring variations to improve NFAc value.
 */
static double ell_ring_improve( PImageDouble angles, PImageInt used, 
                                Ring *ering, double *foci, int x, int y, 
                                int grad_dir, Point **local_buff, 
                                int *size_local_buff, Point **tmp_buff, 
                                int *label, double logNT_ell, FILE *fdebug )
{ 
  double log_nfa, log_nfa_new;
  double delta = 0.5;
  int n;
  int size_tmp_buff;
  Point *tmp;
  fprintf(fdebug,"ell_ring_improve \n");
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

  (*label)++;
  *size_local_buff = 0;
  log_nfa = ell_ring_nfa( angles, used, ering, foci, x, y, grad_dir, 
                          *local_buff, size_local_buff, *label, logNT_ell, fdebug );
     
  /* try to improve the ring, by reducing width; first reduce wmax */
  Ring r; 
  copy_ring( ering, &r );
  for( n=0; n<10; n++ )
    {        
      if( r.wmax-delta > delta )
        {
          r.wmax -= delta;
          size_tmp_buff = 0;
          (*label)++;
          fprintf(fdebug, "wmin %f wmax %f\n",  r.wmin, r.wmax );
          log_nfa_new = ell_ring_nfa( angles, used, &r, foci, x, y, grad_dir,
                                      *tmp_buff, &size_tmp_buff, *label,
                                      logNT_ell, fdebug );
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
  copy_ring( ering, &r );
  for( n=0; n<10; n++ )
    {         
      if( r.wmin + delta < -delta )
        {
          r.wmin += delta;
          size_tmp_buff = 0;
          (*label)++;
          fprintf(fdebug, "wmin %f wmax %f\n",  r.wmin, r.wmax );
          log_nfa_new = ell_ring_nfa( angles, used, &r, foci, x, y, grad_dir,
                                      *tmp_buff, &size_tmp_buff, *label,
                                      logNT_ell, fdebug );
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
/*------------------------------- BUILD OUTPUT -------------------------------*/
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/** A meaningful polygon has been found; store it in output. 
 */
void add_poly_out( Polygon **poly_out, PolyRect *poly, int poly_count, 
                   int *poly_count_max, int **labels, int lbl, FILE *fdebug )
{
  fprintf(fdebug,"add_poly_out poly %d poly_max %d\n", poly_count,*poly_count_max );
  if( poly_count>*poly_count_max )
    {
      *poly_out = (Polygon *) realloc ( *poly_out, (2*poly_count) *
                                                   sizeof(Polygon) );
      *labels = (int *) realloc ( *labels, (2*poly_count) * sizeof(int) );
      *poly_count_max = 2 * poly_count;
    }
  
  polyrect2polygon( poly, &( (*poly_out)[poly_count-1] ) );
  (*labels)[poly_count-1] = lbl;
  fprintf(fdebug,"end add_poly_out \n");
}


/*----------------------------------------------------------------------------*/
/** A meaningful ellipse has been found; store it in output. 
 */
void add_ell_out( Ring **ell_out, Ring *ell, int ell_count, int *ell_count_max,
                  int **labels, int lbl, FILE *fdebug )
{
  fprintf(fdebug,"add_ell_out ell %d ell_max %d \n", ell_count, *ell_count_max);
  if( ell_count>*ell_count_max )
    {
      *ell_out = ( Ring *) realloc ( *ell_out, (2*ell_count) * sizeof(Ring) );
      *labels = ( int *) realloc ( *labels, (2*ell_count) * sizeof(int) );
      *ell_count_max = 2 * ell_count;
    }
  
  copy_ring( ell, &( (*ell_out)[ell_count-1] ) );
  (*labels)[ell_count-1] = lbl;
  fprintf(fdebug,"end add_ell_out \n");
}


/*----------------------------------------------------------------------------*/
/** Entry point in detector's code. Gets as input an image, and returns two 
    lists of geometric shapes: polygons and ellipses (circles are considered 
    particular ellipses, and segments are particular polygons). 
 */
void ELSDc( PImageDouble in, int *ell_count, Ring **ell_out, int **ell_labels,
            int *poly_count, Polygon **poly_out, int **poly_labels, 
            PImageInt out, FILE *fdebug )
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
  char best_type;           /* type of best feature; 'p' for polygon, 
                               'e' for ellipse */
  PImageDouble imgauss;     /* smooth version of the original image, used during
                               candidate selection */
  PImageDouble angles;      /* gradient angles smooth image */
  PImageDouble angles0;     /* gradient angles original image */
  PImageDouble gradmag;     /* gradient magnitude */
  PImageDouble gradx, grady;/* gradient orientations on Ox and Oy */
  PImageInt used;           /* image to mark already used points */
  Rectangle rec;            /* rectangle (line segment) parameters  */
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
  min_size_ell =round( (-logNT_ell-mlog10eps) / log10(p) );

  logNT_seg = 2.5 * (log10((double)xsize) + log10((double)ysize) ) 
                     + 2.0*LOG10_2;
  /* perform gaussian smoothing */
  imgauss = gaussian_sampler( in, 1.0, 0.6 );
  
  /* compute gradient magnitude and orientation  for the smooth image */
  img_gradient_sort( in, rho, &list_p, &mem_p, n_bins, max_grad, &angles, 
                     &gradmag, &gradx, &grady );

  /* compute gradient orientation for the original image */
  angles0 = img_gradient_angle( in, rho );
 
  /* input image not needed any more; free it */
  free_PImageDouble(in);
  
  msize = 100*max(xsize,ysize);

  poly = new_polyrect();
  seg = new_polyrect();

  
  /* begin primitive detection */
  for( ; list_p; list_p = list_p->next )
    { 
      fprintf(fdebug,"loop x %d y %d ecount %d pcount %d \n", 
              list_p->x, list_p->y, *ell_count, *poly_count);
      adr = list_p->y * used->xsize + list_p->x;
      /* only pixels that haven't been used already can be seed pixels */
      if( used->data[adr] == USED ) continue;
      //printf("ok1\n");
      /* only pixels with defined gradient can be seed pixels */
      if( angles->data[adr] == NOTDEF ) continue;
      //printf("ok2\n");
      /* only pixels that are local maxima in the image can be seed pixels */  
      //if( !local_max( imgauss, list_p->x, list_p->y ) ) continue; 
      //printf("ok3\n");
      /* init intra-iteration variables */ 	
      log_nfa = best_nfa = mlog10eps;
      best_feature = NULL; 
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
      fprintf(fdebug,"------------------------------curve_grow start %d \n",reg_size);
      if( !curve_grow( gradmag, angles, used, reg, &reg_size, density_th,
                       prec, poly, &label, pext1, pext2, &spir, fdebug ) ) 
        continue;

      
      /* Validate polygon */
      /* Compute number of tests for polygon; it must be done for each 
         polygon, as it depends on the number of segments in the 
         polygon */
      logNT_poly = ( poly->dim + 1.5 ) * ( log10((double)xsize ) + 
                     log10( (double)ysize) ) + ( poly->dim+1.0 ) * LOG10_2;
      /* polygon's best validation score */
      min_size_poly =(int)((-logNT_poly - mlog10eps)/log10(p));
      fprintf(fdebug,"------------------------------validation start %d \n",reg_size);
      
      if( reg_size > min_size_poly )
        {
          best_nfa = poly_improve( angles0, used, poly, &best_buff, 
                                   &size_best_buff, &tmp_buff, logNT_poly, fdebug );
          best_feature = (void*)poly;
          best_type = 'p';
        }
      fprintf(fdebug,"best poly %d \n",size_best_buff);          

      /* Estimate a single segment on the entire set and validate; 
         this prevents from considering very small regions at the ends 
         of a large segment as new segments. */
      get_seg( angles, gradmag, used, seg, reg, reg_size, prec, fdebug );
      log_nfa = poly_improve( angles0, used, seg, &new_buff, 
                              &size_new_buff, &tmp_buff, logNT_seg, fdebug );
    
      /* if more meaningful than previous, update best feature */
      if( log_nfa > best_nfa ) 
        {
          best_nfa = log_nfa;
          best_feature = (void*)seg;
          best_type = 'p';
          swap_ptr_pts( &best_buff, &new_buff );
          size_best_buff = size_new_buff;
        }
      fprintf(fdebug,"best seg %d \n",size_best_buff);

      /* Estimate circle and ellipse parameters (centre, axes, orientation)
         on the gathered pixels */
      fprintf(fdebug,"reg size before circ valid %d\n",reg_size);
      if( reg_size > min_size_ell ) 
        conic_fit( gradx, grady, reg, reg_size, &buff_fit, &size_buff_fit, 
                   cparam, eparam );

      fprintf(fdebug,"CC\n");
      fprintf(fdebug,"%f %f %f %f %f\n",cparam[0],cparam[1],cparam[2],cparam[3],cparam[4]);
      fprintf(fdebug,"EE\n");
      fprintf(fdebug,"%f %f %f %f %f\n",eparam[0],eparam[1],eparam[2],eparam[3],eparam[4]);
      /* Get parameters of the circle ring and compute validation score */
      ang0 = angles->data[ reg[0].y*angles->xsize + reg[0].x ];
      if( get_ring( reg, reg_size, ang0, cparam, pext1, pext2, CIRCLE,
                    spir, &cring, &grad_dirc, foci ) )
        {
          log_nfa = circ_ring_improve( angles0, used, &cring, reg[0].x, 
                                       reg[0].y, grad_dirc, &new_buff, 
                                       &size_new_buff, &tmp_buff, &label,
                                       logNT_ell, fdebug );
          /* if more meaningful than previous, update best feature */
          if( log_nfa > best_nfa ) 
            {
              best_nfa = log_nfa;
              best_feature = (void*)(&cring);
              best_type = 'c';
              swap_ptr_pts( &best_buff, &new_buff );
              size_best_buff = size_new_buff;
            }
        } /* ! get cring */
      fprintf(fdebug,"best circ %d \n",size_best_buff);   

      /* Get parameters of the ellipse ring and compute validation score */
      if( get_ring( reg, reg_size, ang0, eparam, pext1, pext2, ELLIPSE,
                    spir, &ering, &grad_dire, foci ) )
        {
          log_nfa = ell_ring_improve( angles0, used, &ering, foci, reg[0].x, 
                                      reg[0].y, grad_dire, &new_buff, 
                                      &size_new_buff, &tmp_buff, &label,
                                      logNT_ell, fdebug ); 
          /* if more meaningful than previous, update best feature */
          if( log_nfa > best_nfa ) 
            {
              best_nfa = log_nfa;
              best_feature = (void*)(&ering);
              best_type = 'c';
              swap_ptr_pts( &best_buff, &new_buff );
              size_best_buff = size_new_buff;
            }
        }  /* ! get ering */
      fprintf(fdebug,"best ell %d \n",size_best_buff);

      fprintf(fdebug,"--------------------here type %c nfa %lf\n", best_type, best_nfa);
      /* At this point, we have the most meaninful feature at this iteration. 
         If it is meaningful compared to eps, then store */   
      if( best_nfa <= mlog10eps ) continue;
      fprintf(fdebug, "VALID!!!!!\n");
      label_meaningful++;
      if( best_type == 'c' )
        {
          (*ell_count)++;
          add_ell_out( ell_out, (Ring*)best_feature, *ell_count, &ell_count_max, 
                       ell_labels, label_meaningful, fdebug );
        }
      else
        {
          (*poly_count)++;
          add_poly_out( poly_out, (PolyRect*)best_feature, *poly_count, 
                        &poly_count_max, poly_labels, label_meaningful, fdebug );
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
