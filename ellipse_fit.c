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


  ellipse_fit.c - This file belongs to ELSDc project (Ellipse and Line Segment 
                  Detector with continuous validation).
                - It contains functions to normalize a set of points and to fit 
                  an ellipse or a circle to these points.

------------------------------------------------------------------------------*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <stddef.h>
#include "misc.h"
#include "lapack_wrapper.h"
#include "ellipse_fit.h"


/*----------------------------------------------------------------------------*/
/** Convert ellipse from matrix form to common form:
    ellipse = (centrex,centrey,ax,ay,orientation).
 */
static void ellipse2param( double *p, double *param)
{
  /* p = [ a,     1/2*b, 1/2*d,   
           1/2*b,     c, 1/2*e,  
           1/2*d, 1/2*e, f     ]; */
  double a, b, c, d, e, f;
  double thetarad, cost, sint, cos_squared, sin_squared, cos_sin;
  double Ao, Au, Av, Auu, Avv;
  double tuCentre, tvCentre, wCentre, uCentre, vCentre, Ru, Rv;

  /* check parameters */
  if( p == NULL ) error("ellipse2param: null input ellipse matrix.");
  if( param == NULL ) error("ellipse2param: output 'param' should be non-null");

  a =   p[0];
  b = 2*p[1];
  c =   p[4];
  d = 2*p[2];
  e = 2*p[5];
  f =   p[8]; 

  thetarad = 0.5*atan2(b,a-c); 
  cost = cos(thetarad); 
  sint = sin(thetarad);
  sin_squared = sint * sint;
  cos_squared = cost * cost;
  cos_sin = sint * cost;
  Ao  =  f;
  Au  =  d*cost + e*sint;
  Av  = -d*sint + e*cost;
  Auu = a*cos_squared + c*sin_squared + b*cos_sin;
  Avv = a*sin_squared + c*cos_squared - b*cos_sin;

  if( (Auu == 0) || (Avv==0) )
    { 
      param[0] = 0;
      param[1] = 0;
      param[2] = 0;
      param[3] = 0;
      param[4] = 0;
    }
  else
    {
      tuCentre = -Au / (2.*Auu);
      tvCentre = -Av / (2.*Avv);

      wCentre  =  Ao-Auu*tuCentre*tuCentre - Avv*tvCentre * tvCentre;
      uCentre  =  tuCentre*cost - tvCentre*sint;
      vCentre  =  tuCentre*sint + tvCentre*cost;

      Ru = -wCentre / Auu;
      Rv = -wCentre / Avv;

      if( Ru>0 ) Ru =  pow( Ru,0.5);
      else       Ru = -pow(-Ru,0.5);

      if( Rv>0 ) Rv =  pow( Rv,0.5);
      else       Rv = -pow(-Rv,0.5);

      param[0] = uCentre;
      param[1] = vCentre;
      param[2] = Ru;
      param[3] = Rv;
      param[4] = thetarad;
    }
}


/*----------------------------------------------------------------------------*/
/** Compute normalisation matrix (translates and normalises a set of 
    2D homogeneous points so that their centroid is at the origin and 
    their mean distance from the origin is sqrt(2).
    Input points in 'reg' are in cartesian coordinates (x,y).
    Output matrix 'T' is 3 x 3.
 */
static void vgg_conditioner_from_points( double *T, Point *reg, int reg_size )
{
  double mx = 0.0, my = 0.0;
  double qmean = 0.0, Qt = 0.0, val;
  int i;
  /* check parameters */
  if( reg == NULL ) error("vgg_conditioner_from_points: invalid points list.");
  if( reg_size <= 0 ) error("vgg_conditioner_from_points: invalid size list.");
  
  /* compute mean point */
  for( i=0; i<reg_size; i++ )
    {
      mx += reg[i].x;
      my += reg[i].y;
    }  
  mx /= (double)reg_size; my /= (double)reg_size;

  /* compute mean variance */
  for( i=0; i<reg_size; i++ )
    {
      val = (reg[i].x - mx)*(reg[i].x - mx)+(reg[i].y - my)*(reg[i].y - my);
      Qt += sqrt(val);
    }  	
  qmean = Qt/(double)reg_size;    
  
  /* build normalization matrix */
  val = SQRT2/qmean;
  T[1] = T[3] = 0;
  T[0] = T[4] = val;
  T[2] = -val * mx;
  T[5] = -val * my;
  T[6] = T[7] = 0;
  T[8] = 1;  
} 


/*----------------------------------------------------------------------------*/
/** antisym(u)  A = [ 0,-u(2),u(1); u(2),0,-u(0); -u(1),u(0),0 ]; 
 */
static void antisym( double *u, double *A )
{
  A[0] = A[4] = A[8] = 0;
  A[1] = -u[2];
  A[2] =  u[1];
  A[3] =  u[2];
  A[5] = -u[0];
  A[6] = -u[1];
  A[7] =  u[0];
} 


/*----------------------------------------------------------------------------*/
/** Compute equations for circle/ellipse fitting. Needs buffer memory zone to do 
    the computations.
 */
static void fit_equations( PImageDouble gradx, PImageDouble grady, double *vgg, 
                           Point *reg, int reg_size, double **buff, 
                           int *size_buff_max )
{
  int i,j;
  double K[27];
  double asym[9];
  int idx;
  double crosspr[3]; 
  double pnormx, pnormy, dirnormx, dirnormy;
  int addr;
  /* check parameters */
  if( (gradx == NULL) || (gradx->data == NULL) )
    error("fit_equations: invalid input 'gradx' image.");
    if( (grady == NULL) || (grady->data == NULL) )
    error("fit_equations: invalid input 'grady' image.");
  if( vgg == NULL ) error("fit_equations: invalid normalization matrix.");
  if( reg == NULL ) error("fit_equations: invalid region.");
  if( reg_size<=0 ) error("fit_equations: invalid region size.");
  if( *buff == NULL ) error("fit_equations: invalid buffer.");

  /* if buffer too small, allocate twice the memory required at this step */
  if( reg_size * 24 > *size_buff_max )
    {
      *buff = (double*) realloc( *buff, 
                        sizeof(double) * 2 * reg_size * 24 );  
      if( *buff == NULL ) error("fit_equations: not enough memory.");
      *size_buff_max = 2 * reg_size * 24;
    }


  /* compute normalisation matrix */
  vgg_conditioner_from_points( vgg, reg, reg_size );

  /* compute equation system */
  for( i=0; i<reg_size; i++)
    { 
      idx = i*4*6;
      /* normalise point (pnormx,pnormy) = VGG*(x,y) */
      pnormx = vgg[0]*reg[i].x + vgg[1]*reg[i].y + vgg[2];
      pnormy = vgg[3]*reg[i].x + vgg[4]*reg[i].y + vgg[5];
      
      /* normalise gradient direction (dirnormx,dirnormy) = VGG*(dx,dy) */
      addr = reg[i].y * gradx->xsize + reg[i].x;
      dirnormx = -vgg[0] * grady->data[addr] + vgg[1] * gradx->data[addr];
      dirnormy = -vgg[3] * grady->data[addr] + vgg[4] * gradx->data[addr];

      /* cross product (pnormx,pnormy) x (dirnormx,dirnormy) = tangent line */
      crosspr[0] = -dirnormy;
      crosspr[1] =  dirnormx;
      crosspr[2] =  pnormx * dirnormy - pnormy * dirnormx;

      /* tangent equation: eq = -transpose(kron(TPts(1:3,i),antisym(l)))*J; */
      antisym(crosspr,asym);

      for (j=0;j<9;j++) K[j] = asym[j]*pnormx; 
      for (j=0;j<9;j++) K[j+9] = asym[j]*pnormy; 
      for (j=0;j<9;j++) K[j+18] = asym[j];

      (*buff)[idx]   = -K[0];                      
      (*buff)[idx+1] = -(K[3]+K[9]);      
      (*buff)[idx+2] = -(K[6]+K[18]);       
      (*buff)[idx+3] = -K[12];              
      (*buff)[idx+4] = -(K[15]+K[21]);   
      (*buff)[idx+5] = -K[24];                      

      (*buff)[idx+6]   = -K[1];
      (*buff)[idx+6+1] = -(K[4]+K[10]); 
      (*buff)[idx+6+2] = -(K[7]+K[19]); 
      (*buff)[idx+6+3] = -K[13];
      (*buff)[idx+6+4] = -(K[16]+K[22]);
      (*buff)[idx+6+5] = -K[25];               
      
      (*buff)[idx+12]   = -K[2];      
      (*buff)[idx+12+1] = -(K[5]+K[11]);  
      (*buff)[idx+12+2] = -(K[8]+K[20]);
      (*buff)[idx+12+3] = -K[14];  
      (*buff)[idx+12+4] = -(K[17]+K[23]);  
      (*buff)[idx+12+5] = -K[26];  

      /* position equation: eq = transpose(kron(TPts(1:3,i),TPts(1:3,i)))*J; */
      (*buff)[idx+18] = pnormx * pnormx; 
      (*buff)[idx+19] = 2 * pnormx * pnormy; 
      (*buff)[idx+20] = 2 * pnormx;
      (*buff)[idx+21] = pnormy * pnormy; 
      (*buff)[idx+22] = 2 * pnormy;
      (*buff)[idx+23] = 1;
    }
}


/*----------------------------------------------------------------------------*/
/** Algebraic circle fitting using positional and tangential constraints.
    Equations are already computed and stored in buff.
 */
static void fit_circle( double *buff, double *vgg, int reg_size, double *param )
{
  int i, j, k;
  double A[16]; 
  /* check parameters */
  if( reg_size<=0 ) error("fit_circle: invalid size.");
  if( buff == NULL ) error("fit_circle: invalid buffer.");
  if( vgg == NULL ) error("fit_circle: invalid normalization matrix.");
  if( param == NULL ) error("fit_circle: param must be non null.");

  /* circle fitting computes less parameters, so the equation matrix has 
     only four columns, not six; modify accordingly */
  for( i=0; i<reg_size*4*6; i+=6 )
    { 
      buff[i]   = buff[i] + buff[i+3];
      buff[i+1] = buff[i+2];
      buff[i+2] = buff[i+4];
      buff[i+3] = buff[i+5];
    } 
  /* A = EQ'*EQ; */
  for( i=0; i<16; i++) A[i] = 0.0;

  for( i=0; i<4; i++ )
    for( j=0; j<4; j++)
      for(k=0; k<4*reg_size; k++)
        A[i*4+j] += buff[k*6+i] * buff[k*6+j];

  /* lapack call to solve linear system */
  lap_eig( A, 4 ); 
 
  double s[9];
  s[0] = s[4] = A[0];
  s[1] = s[3] = 0;
  s[2] = s[6] = A[1];
  s[5] = s[7] = A[2];
  s[8] = A[3];

  /* apply inverse(normalisation matrix) */  
  /* C = T'*[ x(1),0,x(2); 0,x(1),x(3) ; x(2),x(3),x(4)]*T;*/

  double C[9];
  C[0] = vgg[0]*vgg[0]*s[0] + vgg[3]*vgg[3]*s[4];
  C[1] = vgg[0]*vgg[1]*s[0] + vgg[3]*vgg[4]*s[4];
  C[2] = vgg[0]*vgg[2]*s[0] + vgg[3]*vgg[5]*s[4]+vgg[0]*s[2] + vgg[3]*s[5];

  C[3] = vgg[0]*vgg[1]*s[0] + vgg[3]*vgg[4]*s[4];
  C[4] = vgg[1]*vgg[1]*s[0] + vgg[4]*vgg[4]*s[4];
  C[5] = vgg[1]*vgg[2]*s[0] + vgg[4]*vgg[5]*s[4] + vgg[1]*s[2] + vgg[4]*s[5];

  C[6] = vgg[0]*vgg[2]*s[0] + vgg[0]*s[6] + vgg[3]*vgg[5]*s[4] + vgg[3]*s[7];
  C[7] = vgg[1]*vgg[2]*s[0] + vgg[1]*s[6] + vgg[4]*vgg[5]*s[4] + vgg[4]*s[7];
  C[8] = vgg[2]*vgg[2]*s[0] + vgg[2]*s[6] + vgg[5]*vgg[5]*s[4] + vgg[5]*s[7]+
         vgg[2]*s[2] + vgg[5]*s[5] + s[8];
  
  ellipse2param( C, param );
}


/*----------------------------------------------------------------------------*/
/** Algebraic ellipse fitting using positional and tangential constraints.
 */
static void fit_ellipse( double *buff, double *vgg, int reg_size, 
                         double *param )
{ 
  int i,j,k;
  double A[36];  
  /* check parameters */
  if( reg_size <= 0 ) error("fit_ellipse: invalid size.");
  if( buff == NULL ) error("fit_ellipse: invalid buffer.");
  if( vgg == NULL ) error("fit_ellipse: invalid normalization matrix.");
  if( param == NULL ) error("fit_ellipse: param must be non null.");
  
  /* A = EQ'*EQ; */
  for( i=0; i<36; i++ ) A[i] = 0.0;

  for( i=0; i<6; i++ )
    for( j=0; j<6; j++ )
      for( k=0; k<4*reg_size; k++ )
        A[i*6+j] += buff[k*6+i] * buff[k*6+j];

  /* lapack call to solve linear system */
  lap_eig( A, 6 ); 
 
  double s[9];
  s[0] = A[0];
  s[1] = s[3] = A[1];
  s[2] = s[6] = A[2];
  s[4] = A[3];
  s[5] = s[7] = A[4];
  s[8] = A[5];
  
  /* apply inverse(normalisation matrix) */
  /* C = T'*[ x(1),x(2),x(3); x(2),x(4),x(5) ; x(3),x(5),x(6)]*T; */
  
  double C[9];
  C[0] = vgg[0]*vgg[0]*s[0] + vgg[0]*vgg[3]*s[3] +
         vgg[0]*vgg[3]*s[1] + vgg[3]*vgg[3]*s[4];
  C[1] = vgg[0]*vgg[1]*s[0] + vgg[1]*vgg[3]*s[3] +
         vgg[0]*vgg[4]*s[1] + vgg[3]*vgg[4]*s[4];
  C[2] = vgg[0]*vgg[2]*s[0] + vgg[2]*vgg[3]*s[3] +
         vgg[0]*vgg[5]*s[1] + vgg[3]*vgg[5]*s[4] + vgg[0]*s[2] + vgg[3]*s[5];

  C[3] = vgg[0]*vgg[1]*s[0] + vgg[0]*vgg[4]*s[3] +
         vgg[1]*vgg[3]*s[1] + vgg[3]*vgg[4]*s[4];
  C[4] = vgg[1]*vgg[1]*s[0] + vgg[1]*vgg[4]*s[3] +
         vgg[1]*vgg[4]*s[1] + vgg[4]*vgg[4]*s[4];
  C[5] = vgg[1]*vgg[2]*s[0] + vgg[2]*vgg[4]*s[3] +
         vgg[1]*vgg[5]*s[1] + vgg[4]*vgg[5]*s[4] + vgg[1]*s[2] + vgg[4]*s[5];

  C[6] = vgg[0]*vgg[2]*s[0] + vgg[0]*vgg[5]*s[3] +
         vgg[0]*s[6] + vgg[2]*vgg[3]*s[1] + vgg[3]*vgg[5]*s[4] + vgg[3]*s[7];
  C[7] = vgg[1]*vgg[2]*s[0] + vgg[1]*vgg[5]*s[3] +
         vgg[1]*s[6] + vgg[2]*vgg[4]*s[1] + vgg[4]*vgg[5]*s[4] + vgg[4]*s[7];
  C[8] = vgg[2]*vgg[2]*s[0] + vgg[2]*vgg[5]*s[3] + vgg[2]*s[6] + 
         vgg[2]*vgg[5]*s[1] + vgg[5]*vgg[5]*s[4] + vgg[5]*s[7] + vgg[2]*s[2] + 
         vgg[5]*s[5] + s[8];
  
  ellipse2param(C,param);
}


/*----------------------------------------------------------------------------*/
/** Given a set of oriented pixels, fit an ellipse and a circle on it
 */
void conic_fit( PImageDouble gradx, PImageDouble grady, Point *reg, 
                int reg_size, double **buff, int *size_buff_max, 
                double *cparam, double *eparam )
{
  double vgg[9];
  if( (gradx == NULL) || (gradx->data == NULL) )
    error("conic_fit: invalid input 'gradx' image.");
  if( (grady == NULL) || (grady->data == NULL) )
    error("conic_fit: invalid input 'grady' image.");
  if( reg == NULL ) error("conic_fit: invalid region.");
  if( reg_size<=0 ) error("conic_fit: invalid region size.");
  if( *buff == NULL ) error("conic_fit: invalid buffer.");
  if( cparam == NULL ) error("conic_fit: invalid 'cparam'.");
  if( eparam == NULL ) error("conic_fit: invalid 'eparam'.");
  
  fit_equations( gradx, grady, vgg, reg, reg_size, buff, size_buff_max );
  fit_ellipse( *buff, vgg, reg_size, eparam );
  fit_circle ( *buff, vgg, reg_size, cparam );
}

