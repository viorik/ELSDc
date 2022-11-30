ELSDc: Ellipse and Line Segment Detector, with Continuous validation
====================================================================

V. Patraucean, P. Gurdjos, R. Grompone von Gioi  
(Corresponding author: Viorica Patraucean vpatrauc@gmail.com)  
Version 0.9, March 2015

ELSDc is an algorithm for joint ellipse and line segment detection in digital
images. It is based on the a-contrario framework and should work without the
need for parameter tuning. The algorithm is described in the following
manuscript:

  "Joint A Contrario Ellipse and Line Detection",
  V. Patraucean, P. Gurdjos, and R. Grompone von Gioi.


Demo
----

An online demo of the algorithm is available, where you can upload your own
images and run ELSDc:

Demo: http://dev.ipol.im/~jirafa/ipol_demo/elsdc/  
(user: demo, pass: demo)


Python interface
----
Available [here](https://github.com/patricknaughton01/ELSDc). Thank you, Patrick!


Requirements
------------

ELSDc requires CLAPACK/CBLAS library for some linear algebra computations.
Version 3.2.1 was used.


Compilation
-----------

You may need to modify the file `Makefile` in the directory `src` to include the
path to liblapack. Then, you can compile by typing `make` in the command
line. This produces the executable called `elsdc`.


Execution
---------

`./elsdc image.pgm` runs ELSDc on the image specified by `image.pgm`. This
version works only with PGM images. This folder contains the image 'shapes.pgm'
for testing purposes.


Output
------

- `output.svg` contains the execution result in SVG format. `shapes_output.svg`
               contains the expected result for the sample image `shapes.pgm`.

- `labels.pgm` after execution, each pixel in this image is labelled with the
               label of the primitive to which it belongs. `shapes_labels.pgm`
               contains the expected result for the sample image `shapes.pgm`.

- `out_ellipse.txt` contains the parameters of the detected circular/elliptical
                   arcs in the form `label x_c y_c a b theta ang_start
                   ang_end`. `shapes_out_ellipse.txt` contains the expected
                   result for the sample image `shapes.pgm`.

- `out_polygon.txt` contains the parameters of the detected line segments,
                   grouped into polygons, defined through contact points, in
                   the form `label number_of_points x1 y1 x2 y2 x3 y3 ...`.
                   `shapes_out_polygon.txt` contains the expected result for
                   the sample image `shapes.pgm`.

The command will print at the end of the execution the numbers of features of
each type found. To check the installation, run `./elsdc shapes.pgm`. The
output should be similar to the one in `shapes_output.svg`, `shapes_labels.pgm`,
`shapes_out_ellipse.txt` and `shapes_out_polygon.txt`. It should find 66
ellipses and 142 polygons. The exact coordinate values may differ slightly due to
rounding error differences in different systems. The execution time for
`shapes.pgm` is about 1.4s on a Dell notebook.


Source Code Files
-----------------

- `main.c`         contains the main() function; entry point into the
                   application, calls IO functions and the detection function.

- `curve_grow.c`   contains functions for region grow (gather neighbour pixels
                   that share the same gradient orientation) and curve grow
                   (gather regions that describe a convex and smooth contour).

- `rectangle.c`    defines a rectangle structure (i.e. segment with width) to
                   approximate the result of region grow.

- `polygon.c`      defines a polygon structure (i.e. collection of rectangles).

- `ring.c`         defines a (circular or elliptical) ring structure.

- `elsdc.c`        contains refinement and validation functions for different
                   types of primitives (ellipse, circle, line segment).

- `ellipse_fit.c`  contains functions to estimate a circle or an ellipse using
                   pixels positions and their gradient orientations.

- `iterator.c`     functions to count the number of aligned pixels inside a
                   rectangle.

- `lapack_wrapper.c` contains wrappers for lapack functions for linear system
                   solve.

- `pgm.c`          IO functions for pgm image format (the only format supported
                   currently).

- `image.c`        defines structures for image representation and functions for
                   gradient computation.

- `gauss.c`        defines a Gaussian kernel and contains functions to performs
                   Gaussian filtering of an image.

- `misc.c`         contains general-purpose functions and constants definitions.

- `svg.c`          functions to write the result in svg format.


Datasets
--------

We make available two datasets together with their ground truth for
quantitative evaluation of ellipse detectors:

**Dataset1_SyntheticCircles** consists of 20 images, 500x500 pixels, containing
non-overlapping and overlapping circles. The images are corrupted with five
different levels of Gaussian noise, with five different noise realisations for
each noise level and for each image, resulting in 500 image instances in total.
The file `coord.txt` contains the ground truth circle parameters in the form:
`pathToimageName.pgm` `xcenter1 ycenter1 radius1 xcenter2 ycenter2 radius2
xcenter3 ycenter3 radius3`. Ten images of pure Gaussian noise are also added
(folder `pureNoise`), where all detections are false positives.

**Dataset2_CalibrationPatterns** contains (in folder `images`) 40 natural
images of calibration patterns that were included in Higuchi et al.'s package
for camera calibration
[link](http://www.ri.cmu.edu/research_project_detail.html?project_id=617&menu_id=261).
The patterns contain coplanar disjoint and concentric circles. The ground truth
was obtained by manually labelling the contours belonging to circles and rings
projections. For each image in folder `images`, there is a corresponding `.txt`
file in folder `ground_truth` with the same name, containing the primitive
parameters in the form: `xcenter ycenter major_axis minor_axis theta`.

**Dataset3_RealImages** contains real images that were used to produce the
results from the paper "Joint A Contrario Ellipse and Line Detection".


License
-------

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU Affero General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option) any
later version.

