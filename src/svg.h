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

  svg.h - This file belongs to ELSDc project (Ellipse and Line Segment 
          Detector with continuous validation).
        - It defines function prototypes to handle (initialize, close) an 
          SVG file, and to write geometric shapes (lines, circles, ellipses, 
          and arcs) in it.         

------------------------------------------------------------------------------*/


#ifndef SVG_H
#define SVG_H

#include "polygon.h"
#include "ring.h"

FILE *init_svg( char *filename, unsigned int xsize, unsigned int ysize );
void fclose_svg( FILE *fsvg );
void write_svg_poly( FILE *fsvg, Polygon *poly );
void write_svg_circ_arc( FILE *fsvg, Ring *cring );
void write_svg_ell_arc( FILE *fsvg, Ring *ering );

#endif

