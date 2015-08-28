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

  svg.c - This file belongs to ELSDc project (Ellipse and Line Segment 
          Detector with continuous validation).
        - It contains functions to handle (initialize, close) an SVG file, and 
          to write geometric shapes (lines, circles, ellipses, and arcs) in it.
          
------------------------------------------------------------------------------*/


#include <stdio.h>
#include <string.h>
#include <math.h>
#include "polygon.h"
#include "ring.h"
#include "svg.h"


/*----------------------------------------------------------------------------*/
/** Open and initialize an SVG file.
 */
FILE *init_svg( char *filename, unsigned int xsize, unsigned int ysize )
{
  FILE *fsvg;

  /* open file */
  fsvg = fopen( filename, "w");
  if( fsvg == NULL ) error("init_svg: unable to open SVG output file.");

  /* write SVG header */
  fprintf( fsvg,"<?xml version=\"1.0\" standalone=\"no\"?>\n");
  fprintf( fsvg,"<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"\n");
  fprintf( fsvg," \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n");
  fprintf( fsvg,"<svg width=\"%upx\" height=\"%upx\" ",xsize,ysize);
  fprintf( fsvg,"version=\"1.1\"\n xmlns=\"http://www.w3.org/2000/svg\" ");
  fprintf( fsvg,"xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n");

  return fsvg;
}


/*----------------------------------------------------------------------------*/
/** Close an SVG file.
 */
void fclose_svg( FILE *fsvg )
{
  /* close SVG file */
  fprintf( fsvg,"</svg>\n");
  if( fclose(fsvg) == EOF )
    error("fclose_svg: unable to close file while writing SVG file.");
}


/*----------------------------------------------------------------------------*/
/** Write polygon in SVG file. For details, see 
    http://www.w3.org/TR/SVG/paths.html#PathDataEllipticalArcCommands 
 */
void write_svg_poly( FILE *fsvg, Polygon *poly )
{
  int i;

  /* write polygon */
  for ( i=0; i<poly->dim; i+=2 )
    {    
      fprintf( fsvg,"<line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" ",
               poly->pts[i].x, poly->pts[i].y, poly->pts[i+1].x,
               poly->pts[i+1].y );
      /* define style: color and width */
      fprintf( fsvg," fill=\"none\" stroke =\"blue\" stroke-width=\"%d\" />\n"
               , 1 );
    }
}


/*----------------------------------------------------------------------------*/
/** Write circular arc in SVG file. For details, see 
    http://www.w3.org/TR/SVG/paths.html#PathDataEllipticalArcCommands 
 */
void write_svg_circ_arc( FILE *fsvg, Ring *cring )
{
  int fa, fs;
  double ang_start, ang_end, C;
  double x1, y1, x2, y2;

  /* consider small arc, in trigonometric sense */
  fa=0; fs=1;
  
  /* compute angles delimiting the arc */
  ang_start = atan2( cring->y1 - cring->cy, cring->x1 - cring->cx); 
  ang_end   = atan2( cring->y2 - cring->cy, cring->x2 - cring->cx);
  
/* if (almost) complete circle, write full circle, not an arc */
  C = M_2__PI * cring->ax;
  if( (cring->full) || ( (angle_diff( ang_start, ang_end ) 
                           < M_2__PI * SQRT2 / C)
                    && (angle_diff_signed( ang_start, ang_end)>0) ) )
    {
      fprintf( fsvg,"<ellipse cx=\"%f\" cy=\"%f\" rx=\"%f\" ry=\"%f\" stroke="
                    "\"red\" stroke-width=\"%d\" fill=\"none\" /> \n",
                    cring->cx, cring->cy, cring->ax, cring->bx, 1 ); 
    }
  else  /* write an arc */ 
    {  
      /* compute starting and ending points of the arc */    
      x1 = cring->ax * cos(ang_start) + cring->cx;
      y1 = cring->ax * sin(ang_start) + cring->cy;
      x2 = cring->ax * cos(ang_end  ) + cring->cx;
      y2 = cring->ax * sin(ang_end  ) + cring->cy;
      if( (double_equal(x1,x2) && double_equal(y1,y2)) || dist(x1,y1,x2,y2) < 2.0 )
        {
          fprintf( fsvg,"<ellipse cx=\"%f\" cy=\"%f\" rx=\"%f\" ry=\"%f\" stroke="
                    "\"red\" stroke-width=\"%d\" fill=\"none\" /> \n",
                    cring->cx, cring->cy, cring->ax, cring->bx, 1 );
          return;
        }
      /* make sure delimiting angles are positive and ordered, to be able to 
         choose between small/big arc */
      if( ang_start < 0 ) ang_start += M_2__PI;
      if( ang_end   < 0 ) ang_end   += M_2__PI;
  
      if( ang_end < ang_start ) ang_end += M_2__PI;

      if( (ang_end - ang_start) > M_PI) fa = 1;
      /* write starting and ending points, axes, orientation, fa, fs */
      fprintf( fsvg,"<path d=\"M %f,%f A%f,%f %f %d,%d %f,%f\"",
               x1, y1, cring->ax, cring->ax, 0.0, fa, fs, x2, y2 );
      /* define style: color and width */
      fprintf(fsvg," fill=\"none\" stroke =\"red\" stroke-width=\"%d\" />\n",1); 
    }
}


/*----------------------------------------------------------------------------*/
/** Write elliptical arc in SVG file. For details, see 
    http://www.w3.org/TR/SVG/paths.html#PathDataEllipticalArcCommands 
 */
void write_svg_ell_arc( FILE *fsvg, Ring *ering )
{
  int fa, fs;
  double ang_start, ang_end;
  double x1, y1, x2, y2;
 
  /* consider small arc, in trigonometric sense */
  fa=0; fs=1;
  
  if( ering->full ) /* if complete ellipse, write full ellipse, not an arc */
    {
    fprintf(fsvg,"<ellipse transform=\"translate(%f %f) rotate(%f)\" rx=\"%f\" "
    "ry=\"%f\" stroke=\"red\" fill=\"none\" stroke-width=\"%d\" />\n",
    ering->cx, ering->cy, (ering->theta)*180/M_PI, ering->ax, ering->bx, 1 );    
    }
  else /* compute limits of the arc and write it */ 
    {  
      /* compute starting point */     
      rosin_point( ering, ering->x1, ering->y1, &x1, &y1 );
      /* compute ending point */
      rosin_point( ering, ering->x2, ering->y2, &x2, &y2 );
      if( (double_equal(x1,x2) && double_equal(y1,y2)) || dist(x1,y1,x2,y2) < 2.0 )
        {
          fprintf(fsvg,"<ellipse transform=\"translate(%f %f) rotate(%f)\" rx=\"%f\" "
          "ry=\"%f\" stroke=\"red\" fill=\"none\" stroke-width=\"%d\" />\n",
          ering->cx, ering->cy, (ering->theta)*180/M_PI, ering->ax, ering->bx, 1 );
          return;
        }
      /* compute angles delimiting the arc */
      ang_start = atan2( y1 - ering->cy, x1 - ering->cx ); 
      ang_end   = atan2( y2 - ering->cy, x2 - ering->cx ); 

      /* make sure delimiting angles are positive and ordered, to be able to 
         choose between small/big arc */
      if( ang_start < 0 ) ang_start += M_2__PI;
      if( ang_end   < 0 ) ang_end   += M_2__PI;
  
      if( ang_end < ang_start ) ang_end += M_2__PI;
      
      if( (ang_end - ang_start ) > M_PI ) fa = 1;
     
      /* write starting and ending points, axes, orientation, fa, fs */
      fprintf( fsvg,"<path d=\"M %f,%f A%f,%f %f %d,%d %f,%f\"",
               x1, y1, ering->ax, ering->bx, ering->theta*180/M_PI, 
               fa, fs, x2, y2);
      /* define style: color and width */
      fprintf(fsvg," fill=\"none\" stroke =\"red\" stroke-width=\"%d\" />\n",1);   
    }
}
