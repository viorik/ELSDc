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
#include "rectangle.h"
#include "iterator.h"
#include "polygon.h"
#include "ring.h"
#include "curve_grow.h"



/*----------------------------------------------------------------------------*/
/** Group in 'reg' neighbour pixels sharing the same orientation up to 
    precision 'prec'.
    reg_index: 0.....start.....idx_buf (end)
    (start...idx_buf) is the buffer zone (>=1) of the list; the orientations of 
    the pixels in this zone are used to initialize reg_angle. 
 */
static void region_grow( PImageDouble angles, PImageInt used, Point *reg, 
                         int start, int idx_buff, int *end, int label, 
                         int label_start, double prec, double *reg_angle )
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
  /* if( start >= *end ) error("region_grow: invalid indexes.");*/
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
                if( (ang != NOTDEF ) && (is_aligned( ang, *reg_angle, prec )) )
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
                         double cx, double cy, double reg_angle, double prec )
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
void region2rect( PImageDouble gradmag, Point *reg, int start, int end,
                         double reg_angle, Rectangle *rec, double prec )
{
  double cx, cy, dx, dy, l, w, theta;
  double weight,sum,l_min,l_max,w_min,w_max;
  int i, adr;
  /* Check parameters */
  if( (gradmag == NULL) || (gradmag->data == NULL) )
    error("region2rect: invalid image 'gradmag'.");
  if( reg == NULL ) error("region2rect: invalid input region.");
  if( start >= end ) error("region2rect: invalid indexes."); 
  if( rec == NULL ) error("region2rect: invalid 'rec'.");

  /* Center */
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

  /* Orientation */  
  theta = get_theta( gradmag, reg, start, end, cx, cy, reg_angle, prec );

  /* Length and width */
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

  /* Store values */
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
  rec->len = dist( rec->x1, rec->y1, rec->x2, rec->y2 );
  /* Correct if rectangle too thin */
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
  /* Check parameters */
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

  /* Compute region points density */
  density = (double) (*end - start) /
                     (dist(rec->x1,rec->y1,rec->x2,rec->y2) * rec->width);

  if( density >= density_th ) return TRUE;

  /* Compute region radius */
  rad1 = dist( (double)x, (double)y, rec->x1, rec->y1 );
  rad2 = dist( (double)x, (double)y, rec->x2, rec->y2 );
  rad = rad1 > rad2 ? rad1 : rad2;

  while( density < density_th )
    {
      rad *= 0.75;

      /* Remove points from the region and update 'used' map */
      for( i=start; i<*end; i++ )
        if( dist( x, y, (double) reg[i].x, (double) reg[i].y ) > rad )
          {
            /* Point not kept, mark it as NOTUSED */
            used->data[ reg[i].x + reg[i].y * used->xsize ] = NOTUSED;
            /* Remove point from the region */
            reg[i].x = reg[*end-1].x; /* if i==*end-1 copy itself */
            reg[i].y = reg[*end-1].y;
            --(*end);
            --i; /* to avoid skipping one point */
          }
      /* Reject if the region is too small.
         2 is the minimal region size for 'region2rect' to work. */

      if( (*end - start) < 2 ) return FALSE;

      /* Re-compute rectangle */
      region2rect( gradmag, reg, start, *end, reg_angle, rec, prec );

      /* Re-compute region points density */
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
                   double prec )
{
  double angle, ang_d, mean_angle, tau, density, ang_c, sum, s_sum, reg_angle;
  double gradmax;
  int i, n, xx, yy;
  int xc, yc;
  
  /* Check parameters */
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
 
  /* Compute region points density */
  density = (double) (*end - start) /
            ( dist(rec->x1,rec->y1,rec->x2,rec->y2) * rec->width );

  if( density >= density_th ) return TRUE;
  
  /*------ First try: reduce angle tolerance ------*/

  /* Compute the new mean angle and tolerance */
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

 
  /* Find a new region from the same starting point and new angle tolerance */
  *end = idx_buff;
  (*label)++;
  region_grow( angles, used, reg, start, idx_buff, end, *label, label_start, 
               tau, &reg_angle );

  /* If the region is too small, reject */
  if( *end - start <= 2 ) return FALSE;

  /* Re-compute rectangle */
  region2rect( gradmag, reg, start, *end, reg_angle, rec, prec );

  /* Re-compute region points density */
  density = (double) (*end - start) /
                     ( dist(rec->x1,rec->y1,rec->x2,rec->y2) * rec->width);

  /*------ Second try: reduce region radius ------*/
  if( density < density_th )
    return reduce_region_radius( gradmag, angles, used, reg, start, end, 
                                 reg_angle, rec, xc, yc, density_th, prec );
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
                     double *pext, double prec, int label, int label_start )
{
  int xx, yy, i;
  int adr;
  int idx_buff0 = *idx_buff; 
  double width;
  double l1, l2, l;
  
  /* Check parameters */
  if( (angles == NULL) || (angles->data == NULL) )
    error("px_seed: invalid input 'angles' image.");
  if( (used == NULL) || (used->data == NULL) )
    error("px_seed: invalid input 'used' image.");
  if( rec == NULL ) error("px_seed: invalid input rectangle.");
  if( reg == NULL ) error("px_seed: invalid input region.");
  if( end <= start ) error("px_seed: invalid region size.");

  width = rec->width/2.0;

  /* Distance between rectangle centre and extreme point; a new point can be 
     a seed only if it is at a distance farther than this */ 
  l1 = fabs( (pext[0] - rec->cx) * rec->dx + (pext[1] - rec->cy) * rec->dy);

  /* Scan previous region to find valid points that are close to the end of 
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
              /* If point in image */
              if( in_image( xx, yy, used->xsize, used->ysize ) )
                {                   
                  adr = yy * used->xsize + xx;                  
                  /* If point not already used or visited or notdef */
                  if( (used->data[adr] != USED) && (used->data[adr] != label ) &&
                      (angles->data[adr] != NOTDEF) && used->data[adr] < label_start )
                    {
                      l2 = fabs( (xx - rec->cx) * rec->dx + 
                                 (yy - rec->cy) * rec->dy ); 
                      /* If point's orientation is similar and the point is at
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

  /* If none of the visited pixels satisfy the conditions, take as seed 
     the end of the rectangle, if it is inside the image */
  if( idx_buff0 == *idx_buff )
    {
      xx = floor(pext[0]); yy = floor(pext[1]);
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
                      int *end, double *pext, double *spir, double density_th,
                      double prec,int *sgn, int *label, int label_start, int w )
{
  int convex = 1;
  double reg_angle, difang, ang_prev;
  double lenrec, lenth;
  int buf_tmp = 0;
  Rectangle rec;
  
  /* Check parameters */
  if( (angles == NULL) || (angles->data == NULL) )
    error("subcurve: invalid input 'angles' image.");
  if( (gradmag == NULL) || (gradmag->data == NULL) )
    error("subcurve: invalid input 'gradmag' image.");
  if( (used == NULL) || (used->data == NULL) )
    error("subcurve: invalid input 'used' image.");
  if( reg == NULL ) error("subcurve: invalid input region.");
  if( *end <= start ) error("subcurve: invalid region size.");
  if( *label <= 1 ) error("subcurve: forbidden label value.");


  /* Get length and angle of initial rectangle (stored in polygon on first 
     position) */
  lenth = poly->rectlist[0].len;
  ang_prev = poly->rectlist[0].theta;
 
  /* Add new rectangles while the conditions of smoothness and convexity 
     are satisfied */
  while( convex )
    {
      (*label)++;
      start = *end;
      *end = idx_buff;
      buf_tmp = *end - start;
      /* Start region grow at the end of the previous rectangle */
      region_grow( angles, used, reg, start, idx_buff, end, *label, label_start,
                   prec, &reg_angle );       
      if( *end - start > buf_tmp )
        {
          /* Estimate and refine rectangle */
          region2rect( gradmag, reg, start, *end, reg_angle, &rec, prec );
          (*label)++;
          if( !refine( gradmag, angles, used, reg, start, idx_buff, end,
                       label, label_start, &rec, density_th, prec) ) 
            {
              /* If operation not successful, clean last visited points and 
                 return */
              mark_img_pts( used, reg, start, *end, NOTUSED );
              *end = start;
              break;
            }
          if( *end - start <= buf_tmp ) 
            {
              *end = start;
              break;
            }

          /* Convexity rule imposes that consecutive rectangles turn in the same
             direction, hence must have the same sign for angle difference 
             between consecutive rectangles' orientation. If set, verify sign
             before adding a new rectangle. If not set (at first call), 
             then set sign. */
          if( *sgn == 0 ) 
            *sgn = sign( angle_diff_signed( rec.theta, ang_prev ) ); 

          /* Get angle difference between current rectangle and previous 
             rectangle */  
          difang = angle_diff_signed( rec.theta, ang_prev );

          /* If consecutive rectangles don't have comparable lengths  
             (condition for smoothness): clean and return; 
             Special case: last rectangle that joins the ends of a full 
             chain can be shorter than previous.  */
          lenrec = rec.len;
          if( lenrec > lenth*3.0 || 
              (lenrec < lenth/3.0 && (*spir + fabs(difang)) <= 0.9*M_2__PI) ) 
            {
              mark_img_pts( used, reg, start, *end, NOTUSED );
              *end = start;
              break;
            }
          lenth = lenrec;


          /* Convexity & non-spiral check */
          if( (sign(difang) == *sgn) && (fabs(difang) < (M_PI/3.0) ) &&
              ( (*spir = *spir + fabs(difang)) <= M_2__PI) )
            { 
              /* If all conditions are met, find next seed point, 
                 and add rectangle to polygon */
              (*label)++;
              idx_buff = *end;
              if( w==0 ) {pext[0] = rec.x1; pext[1] = rec.y1;}
              else {pext[0] = rec.x2; pext[1] = rec.y2;}
              px_seed( angles, used, &rec, reg, start, &idx_buff, *end,
                       pext, prec, *label, label_start );

              /* Add rectangle to polygon */
              poly->dim++;
              add_rect_to_polyrect( poly, &rec );
              ang_prev = rec.theta; 
            }      
          else /* Not convex */
            {
              convex = 0;
              /* Subtract the contribution of the last rectangle and clean */
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
int curve_grow( PImageDouble gradmag, PImageDouble angles, 
                       PImageInt used, Point *reg, int *reg_size, 
                       double density_th, double prec, PolyRect *poly, 
                       int *label, double *pext1, double *pext2, double *spir ) 
{
  int sgn = 0;
  Rectangle rec;
  int size_first_reg = 0;
  int label_start;
  int start, idx_buff;
  double reg_angle;
  
  /* Check parameters */
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
  
  /* Init indexes */
  start = 0;
  idx_buff = *reg_size;
  label_start = (*label);
  /* Perform initial region growing, to estimate first rectangle */
  region_grow( angles, used, reg, start, idx_buff, reg_size, *label, label_start,
               prec, &reg_angle ); 

  if( *reg_size <= 2 ) return 0;
  
  /* Estimate initial rectangle */
  region2rect( gradmag, reg, start, *reg_size, reg_angle, &rec, prec );
  
  /* Refine rectangle if density of aligned points is less than density 
     threshold */
  if( !refine( gradmag, angles, used, reg, start, idx_buff, reg_size, label, 
               label_start, &rec, density_th, prec ) )
    {
      mark_img_pts( used, reg, start, *reg_size, NOTUSED );
      return 0;
    }
  if( rec.len <= 2 ) return 0;

  /* Save size of first region; useful to compute second seed, after scanning
     first direction */
  size_first_reg = *reg_size;

  /* Add rectangle to poly */
  poly->dim++; 
  add_rect_to_polyrect( poly, &rec );

  /* Compute first seed with the first end of the initial rectangle */
  pext1[0] = rec.x2; 
  pext1[1] = rec.y2;
  (*label)++;
  idx_buff = *reg_size;
  px_seed( angles, used, &rec, reg, start, &idx_buff, *reg_size, pext1, prec,
           *label, label_start ); 

  /* Lauch scan in first direction */ 
  subcurve( angles, gradmag, used, poly, reg, start, idx_buff, reg_size, 
            pext1, spir, density_th, prec, &sgn, label, label_start, 1 ); 
  
  /* Scan the other direction only if a complete tour was not scanned; this
     condition is to avoid spiral curves */
  pext2[0] = rec.x1; 
  pext2[1] = rec.y1; 
  if( *spir < M_2__PI )
    {
      /* Compute second seed with the second end of the initial rectangle */  
 
      idx_buff = *reg_size;
      (*label)++;
      px_seed( angles, used, &rec, reg, 0, &idx_buff, size_first_reg, 
               pext2, prec, *label, label_start );
      /* When scanning the other direction, the sign must be inverted */ 
      sgn = -sgn;
      /* Lauch scan in second direction */ 
      subcurve( angles, gradmag, used, poly, reg, start, idx_buff, reg_size, 
                pext2, spir, density_th, prec, &sgn, label, label_start, 0 );
    }    
  return 1;
}
