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


  main.c - This file belongs to ELSDc project (Ellipse and Line Segment 
           Detector with continuous validation)
         - It contains the main function, which reads a pgm image, calls the 
           detector, and writes the result in ASCII and SVG form.

------------------------------------------------------------------------------*/


#include <stdio.h>
#include <stdlib.h>
#include "misc.h"
#include "pgm.h"
#include "svg.h"
#include "polygon.h"
#include "ring.h"
#include "elsdc.h"

/*----------------------------------------------------------------------------*/
/** Main function; takes as argument the name of a pgm image, runs the detection
    on it, and writes the result in ASCII and SVG form.  
 */
int main( int argc, char **argv )
{
  /* check arguments */
  if( argc < 2 ) error("usage: ./elsdc image_name.pgm");

  PImageDouble in;     /* input image */
  PImageInt    out;    /* output image having the same size as 'in'; the pixels
                          supporting a certain geometric primitive are marked 
                          with the same label */

  int ell_count = 0;   /* number of detected ellipses */
  int *ell_labels=NULL;/* the pixels supporting a certain ellipse are marked 
                          with the same unique label */
  Ring *ell_out = NULL;/* array containing the parameters of the detected 
                          ellipses; correlated with ell_labels, i.e. the i-th
                          element of ell_labels is the label of the pixels 
                          supporting the ellipse defined by the parameters 
                          ell[i] */
                       
  int poly_count = 0;  /* number of detected polygons */
  int *poly_labels=NULL;/* the pixels supporting a certain polygon are marked 
                          with the same unique label */
  Polygon *poly_out=NULL;/* array containing the parameters of the detected 
                          polygons; correlated with poly_labels, i.e. the i-th
                          element of ell_labels is the label of the pixels 
                          supporting the polygon defined by the parameters
                          poly[i] */

  FILE *ell_ascii;     /* output file with the parameters of the detected 
                          ellipses -- ASCII format */
  FILE *poly_ascii;    /* output file with the parameters of the detected
                          polygons -- ASCII format */
  FILE *fsvg;          /* output file with the detected ellipses and polygons 
                          in vectorial form */
  int i,j;

  /* read input image; must be PGM form */
  in = read_pgm_image_double( argv[1] );
  int xsize = in->xsize, ysize = in->ysize;
  
  /* create and initialize with 0 output label image */
  out = new_PImageInt_ini( in->xsize, in->ysize, 0 );  
  
  /* call detection procedure */
  ELSDc( in, &ell_count, &ell_out, &ell_labels, &poly_count, &poly_out, 
         &poly_labels, out );

  /* write results in ASCII */
  /* Ellipse file: each line contains 1 integer and 11 doubles in the form
     ell_label x1 y1 x2 y2 cx cy ax bx theta ang_start ang_end
     where:
     ell_label          -- the label of the pixels supporting this ellipse in 
                           'out' labels image 
     x1, y1, x2, y2     -- points delimiting the arc in trigonometric sense
     cx, cy             -- ellipse centre
     ax, bx             -- ellipse axes; ax is the greater axis
     theta              -- orientation of the ellipse
     ang_start, ang_end -- angles delimiting the arc in trigonometric sense */
  if( ell_out != NULL )
    {
      if( (ell_ascii = fopen("out_ellipse.txt","w")) == NULL )
        error("main: can't open ellipse output file.");
      for( i=0; i<ell_count; i++ )
        {
          fprintf( ell_ascii,"%d ", ell_labels[i] );
          fprintf( ell_ascii,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n",
                   ell_out[i].x1, ell_out[i].y1, ell_out[i].x2, ell_out[i].y2, 
                   ell_out[i].cx, ell_out[i].cy, ell_out[i].ax, ell_out[i].bx,
                   ell_out[i].theta, ell_out[i].ang_start, ell_out[i].ang_end );
        }
      fclose(ell_ascii);
    }    

  /* Polygon file: each line contains 2 integers and a variable number of  
     doubles in the form 
     poly_label n x1 y1 x2 y2 ... xn yn
     where:
     poly_label      -- the label of the pixels supporting this polygon in 
                        'out' labels image
     n               -- the number of points in the polygon; it is the double 
                        of the number of segments
     x1 y1 ... xn yn -- (x,y) coordinates of the ending points of each segment.
     A polygon with n/2 line segments has n points, given in consecutive order.
   */
  if( poly_out != NULL )
    {
      if( (poly_ascii = fopen("out_polygon.txt","w")) == NULL )
        error("main: can't open polygon output file.");
      for( i=0; i<poly_count; i++ )
        {
          fprintf( poly_ascii, "%d %d ", poly_labels[i], poly_out[i].dim );
          for( j=0; j<poly_out[i].dim; j++ )
            fprintf( poly_ascii,"%lf %lf ", poly_out[i].pts[j].x, poly_out[i].pts[j].y);
          fprintf( poly_ascii, "\n" );
        }
      fclose(poly_ascii);
    }  

  /* write vectorial output in SVG format */
  if( (ell_out != NULL) || (poly_out != NULL) )
    {
      /* init svg file */
      fsvg = init_svg( "output.svg", xsize, ysize );
 
      /* write ellipses */
      for( i=0; i<ell_count; i++)
        /* distinguish between circle and ellipse, because the procedures 
           to write them are different */
        if( (double_equal( ell_out[i].ax, ell_out[i].bx )) )/* && (ell_out[i].theta == 0) ) */
          write_svg_circ_arc( fsvg, &ell_out[i] );
        else 
          write_svg_ell_arc( fsvg, &ell_out[i] );
        
      /* write polygons */
      for( i=0; i<poly_count; i++ )
        write_svg_poly( fsvg, &poly_out[i] );
      /* close svg file */
      fclose_svg( fsvg );
    }

  /* write labels image in pgm form */
  write_pgm_image_int( out->data, out->xsize, out->ysize, "labels.pgm" );
  free_PImageInt(out);
  if( ell_out != NULL ) {free(ell_out); free(ell_labels);}
  if( poly_out != NULL ) 
    {
      for( i=0; i<poly_count; i++ )
        free(poly_out[i].pts);
      free(poly_out); 
      free(poly_labels);
    }
  printf("Number of ellipses detected = %d\n",ell_count);
  printf("Number of polygons detected = %d\n",poly_count);
  return 0;
}
