<b>ELSDc -- Ellipse and Line Segment Detector, with Continuous validation</b>

The source files included in this repository (folder 'src') contain the source 
code of ELSDc (Ellipse and Line Segment Detector, with Continuous validation),
described in 'Joint A Contrario Ellipse and Line Detection', by 
V. Patraucean, P. Gurdjos, R. Grompone von Gioi; this is the enhanced version 
of ELSD, published in 
'A Parameterless Line Segment and Elliptical Arc Detector with Enhanced 
Ellipse Fitting', V. Patraucean, P. Gurdjos, R. Grompone von Gioi, ECCV2012.

<b>Corresponding author</b>: viorica patraucean vpatrauc@gmail.com

Test online the detector by uploading your own images at 
http://dev.ipol.im/~jirafa/ipol_demo/elsdc/ (user: demo, pass:demo).

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU Affero General Public License as published by the Free 
Software Foundation, either version 3 of the License, or (at your option) any
later version. 

<b>REQUIREMENTS</b>

ELSDc requires CLAPACK/CBLAS library for some linear algebra computations. 
Version 3.2.1 was used.

<b>COMPILATION</b>

Modify in 'makefile' the path to liblapack if needed, and compile by typing <b>make</b>
in the command line. This produces the executable called 'elsdc'.


<b>EXECUTION</b>
<b>./elsdc imagename</b>: runs ELSDc on the image specified by 'imagename'. This 
                   version works only with PGM images. This folder contains the
                   image 'shapes.pgm' for testing purposes.


<b>OUTPUT</b>

<b>output.svg:</b>       contains the execution result in SVG format. 'shapes_output.svg' 
                   contains the result for the sample image 'shapes.pgm'.
<b>labels.pgm:</b>       after execution, each pixel in this image is labelled with the 
                   label of the primitive to which it belongs.   
<b>out_ellipse.txt:</b>  contains the parameters of the detected circular/elliptical 
                   arcs in the form 'label x_c y_c a b theta ang_start ang_end'.
<b>out_polygon.txt:</b>  contains the parameters of the detected line segments, grouped 
                   into polygons, defined through contact points, in the form 
                   'label number_of_points x1 y1 x2 y2 x3 y3 ...'. 

In the console, the numbers of features of each type are displayed. 
To check the installation, run ./elsdc for 'shapes.pgm' image. The output should be
similar to 'shapes_output.svg', and contain 66 ellipses and 145 polygons, whose 
parameters are contained in 'out_ellipse.txt' and 'out_polygon.txt'.
The execution time for this sample image is about 4s on my Dell notebook.  


<b>SOURCE CODE FILES (folder 'src')</b>

<b>main.c:</b>           contains the main() function; entry point into the application,
		   calls IO functions and the detection function.
<b>curve_grow.c:</b>	   contains functions for region grow (gather neighbour pixels 
                   that share the same gradient orientation) and curve grow (gather
                   regions that describe a convex and smooth contour).
<b>rectangle.c:</b>	   defines a rectangle structure (i.e. segment with width) to 
		   approximate the result of region grow.
<b>polygon.c:</b>  	   defines a polygon structure (i.e. collection of rectangles).
<b>ring.c:</b>	   defines a (circular or elliptical) ring structure. 
<b>elsdc.c:</b>	   contains refinement and validation functions for different 
                   types of primitives (ellipse, circle, line segment). 
<b>ellipse_fit.c:</b>	   contains functions to estimate a circle or an ellipse using 
                   pixels positions and their gradient orienatations.
<b>iterator.c:</b>	   functions to count the number of aligned pixels inside a rectangle.
<b>lapack_wrapper.c:</b> contains wrappers for lapack functions for linear system solve.
<b>pgm.c:</b>            IO functions for pgm image format (the only format supported 
                   currently).
<b>image.c:</b>	   defines structures for image representation and functions for 
                   gradient computation.
<b>gauss.c:</b>	   defines a Gaussian kernel and contains functions to performs 
                   Gaussian filtering of an image.
<b>misc.c:</b>	   contains general-purpose functions and constants definitions.
<b>svg.c:</b>     	   functions to write the result in svg format.

<b>DATASETS</b>

We make available two datasets together with their ground truth for 
quantitative evaluation of ellipse detectors:

<b>Dataset1_SyntheticCircles</b> consists of 20 images, 500x500 pixels, containing 
non-overlapping and overlapping circles. The images are corrupted with five 
different levels of Gaussian noise, with five different noise realisations for 
each noise level and for each image, resulting in 500 image instances in total. 
The file 'coord.txt' contains the ground truth circle parameters in the form:
'pathToimageName.pgm'
'xcenter1 ycenter1 radius1 xcenter2 ycenter2 radius2 xcenter3 ycenter3 radius3'.
Ten images of pure Gaussian noise are also added (folder 'pureNoise/'), where 
all detections are false positives.

<b>Dataset2_CalibrationPatterns</b> contains (in folder 'images') 40 natural images 
of calibration patterns that were included in Higuchi et al.’s package for 
camera calibration 
http://www.ri.cmu.edu/research_project_detail.html?project_id=617&menu_id=261. 
The patterns contain coplanar disjoint and concentric circles. The ground truth 
was obtained by manually labelling the contours belonging to circles and rings
projections. For each image in folder 'images', there is a correponding .txt 
file in folder 'ground_truth' with the same name, containing the primitive 
parameters in the form:
xcenter ycenter major_axis minor_axis theta.

<b>Dataset3_RealImages</b> contains real images that were used to produce the results 
from the paper "A Contrario Joint Ellipse and Line Detection".

