/*------------------------------------------------------------------------------

  Copyright (c) 2014 viorica patraucean (vpatrauc@gmail.com)
  
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


  lapack_wrapper.c - This file belongs to ELSDc project (Ellipse and Line 
                     Segment Detector with continuous validation).
                   - It contains wrappers for lapack functions.

------------------------------------------------------------------------------*/

#include <stddef.h>
#include "lapack_wrapper.h"

typedef long int integer;
typedef double doublereal;


/*----------------------------------------------------------------------------*/
/** Header of lapack function for computing the eigen-decomposition of a
    matrix of double. 
 */
int dsyev_(char *jobz, char *uplo, integer *n, doublereal *a, 
           integer *lda, doublereal *w, doublereal *work, integer *lwork, 
           integer *info);


/*----------------------------------------------------------------------------*/
/** Solve linear system. 
    'n' number of unknowns;
    'A' coefficients of the equations. 
 */
void lap_eig( double *A, int n ) 
{
  char jobz = 'V';
  char uplo = 'U';
  integer M = (integer)n;
  integer LDA = M;
  integer LWORK = 24;
  integer INFO;
  doublereal W[6];
  doublereal WORK[LWORK];

  /* Solve eigenproblem */
  dsyev_( &jobz, &uplo, &M, (doublereal*)A, &LDA, W, WORK, &LWORK, &INFO ); 
}

