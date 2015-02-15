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

  pgm.c - This file belongs to ELSDc project (Ellipse and Line Segment 
          Detector with continuous validation).
        - It contains functions to read and write pgm images, and static 
          functions called by these to read characters and numbers from a file.

------------------------------------------------------------------------------*/


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>
#include <ctype.h>
#include <time.h>
#include <limits.h>
#include "misc.h"
#include "image.h"
#include "pgm.h"


/*----------------------------------------------------------------------------*/
/** Skip white characters and comments in a PGM file.
 */
static void skip_whites_and_comments( FILE *f )
{
  int c;
  do
    {
      while( isspace( c = getc(f) ) ); /* skip spaces */
      if( c == '#') /* skip comments */
        while( (c != '\n') && (c != '\r') && (c != EOF) )
          c = getc(f);
    }

  while( (c == '#') || (isspace(c)) );
  if( (c != EOF) && (ungetc( c, f ) == EOF) )
    error("skip_whites_and_comments: unable to 'ungetc' in reading PGM file.");
}


/*----------------------------------------------------------------------------*/
/** Read an ASCII number from a PGM file.
 */
static unsigned int get_num( FILE *f )
{
  unsigned int num;
  int c;

  while( isspace( c = getc(f) ) ); /* skip spaces */
  if( !isdigit(c) ) error("Error: corrupted PGM file.");

  num = (unsigned int) (c - '0');
  while( isdigit( c = getc(f) ) ) num = 10 * num + c - '0';
  if( (c != EOF) && (ungetc( c, f ) == EOF) )
    error("get_num: unable to 'ungetc' while reading PGM file.");

  return num;
}


/*----------------------------------------------------------------------------*/
/** Read a PGM file into an "image_double".
    If the name is "-" the file is read from standard input.
 */
PImageDouble read_pgm_image_double( char *name )
{
  FILE *f;
  int c, bin = FALSE;
  unsigned int xsize, ysize, depth, x, y;
  PImageDouble image;

  /* open file */
  f = fopen( name,"rb");
  if( f == NULL ) 
    error("read_pgm_image_double: unable to open input image file.");

  /* read header */
  if( getc(f) != 'P') error("read_pgm_image_double: not a PGM file.");
  if( (c = getc(f) ) == '2') bin = FALSE;
  else if( c == '5') bin = TRUE;
  else error("read_pgm_image_double: not a PGM file.");

  skip_whites_and_comments(f);
  xsize = get_num(f);            /* X size */
  skip_whites_and_comments(f);
  ysize = get_num(f);            /* Y size */
  skip_whites_and_comments(f);
  depth = get_num(f);            /* depth */

  if( depth==0 ) 
    fprintf( stderr,"Warning: depth=0, probably invalid PGM file.\n");

  /* white before data */
  if( !isspace( c = getc(f) ) ) 
    error("read_pgm_image_double: corrupted PGM file.");

  /* get memory */
  image = new_PImageDouble( xsize, ysize );
  
  /* read data */
  for( y=0; y<ysize; y++)
    for( x=0; x<xsize; x++)
      image->data[ x + y * xsize ] = bin ? (double) getc(f)
                                         : (double) get_num(f);

  /* close file if needed */
  if( (f != stdin) && (fclose(f) == EOF) )
    error("read_pgm_image_double: unable to close file.");
  
  return image;
}


/*----------------------------------------------------------------------------*/
/** Write an int image into a PGM file.
    If the name is "-" the file is written to standard output.
 */
void write_pgm_image_int( int *image, int xsize, int ysize, char *name )
{
  FILE *f;
  int x, y, n, v, max, min;

  /* check input */
  if( (image == NULL) || (xsize <= 0) || (ysize <= 0) )
    error("write_pgm_image_int: invalid input image.");

  /* check min and max values */
  max = min = 0;
  for( y=0; y<ysize; y++ )
    for( x=0; x<xsize; x++ )
      {
        v = image[ x + y * xsize ];
        if( v > max ) max = v;
        if( v < min ) min = v;
      }

  if( min < 0 ) fprintf( stderr,
    "Warning: write_pgm_image_int: negative values in '%s'.\n", name );
  if( max > 65535 ) fprintf( stderr,
    "Warning: write_pgm_image_int: values exceeding 65535 in '%s'.\n", name );

  /* open file */
  if( strcmp( name,"-") == 0 ) f = stdout;
  else f = fopen( name,"w");
  if( f == NULL ) 
    error("write_pgm_image_int: unable to open output image file.");

  /* write header */
  fprintf( f,"P2\n");
  fprintf( f,"%d %d\n", xsize, ysize );
  fprintf( f,"%d\n", max);

  /* write data */
  for( n=0, y=0; y<ysize; y++ )
    for( x=0; x<xsize; x++ )
      {
        fprintf( f,"%d ", image[ x + y * xsize ] );
        if( ++n == 8 )  /* lines should not be longer than 70 characters */
          {
            fprintf( f,"\n");
            n = 0;
          }
      }

  /* close file if needed */
  if( (f != stdout) && (fclose(f) == EOF) )
    error("write_pgm_image_int: unable to close file.");
}
