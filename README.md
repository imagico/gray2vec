
# gray2vec grayscale image vectorizer using pixel aligned halftoning

This tool produces a polygon set from a grayscale image using halftoning.

The polygons are pixel aligned so when rendered with AGG based renderers 
(as well as many other rendering engines) the result will approximately match 
the input image being used as a coverage mask.

The output polygon grid has half the linear resolution of the input image.

The process used tries to minimize data size while maintaining consistency 
with the input data in terms of pixel coverage.  It also tries to maintain 
some similarity regarding the topological configuration (as present in the double
resolution input data).

![Grayscale source raster](http://www.imagico.de/map/pict/sample_fractions.png)

![Generated vector data](http://www.imagico.de/map/pict/sample_vectorized.png)


## Compiling the program

Internally the tool uses some code and components from GDAL 2.  To allow 
building it with GDAL 1 the source file `gdalrasterpolygonenumerator.cpp`
(from GDAL 2.1.1) is included in the repository - building and linking this
should only be required for GDAL 1.x

Dependecies: [GDAL](http://gdal.org/) and [CImg](http://cimg.eu/).


## Program options

The program offers the following command line options:

* `-i` input image (required)
* `-o` output vector file (required)
* `-c` combined input image - to generate accurate output with several layers (optional)
* `-l` output layer name. Default: `polygons`.
* `-x` x attribute to apply to generated polygons (integer value, optional)
* `-y` y attribute to apply to generated polygons (integer value, optional)
* `-z` z attribute to apply to generated polygons (integer value, optional)
* `-complement` process complement (inverse) of input.  Default: `off`.
* `-append` append to existing output file.  Default: `off`.
* `-me` maximum error to accept for pixel coverage fraction.  Default: `0.05`.
* `-debug` generate additional debug output.  Default: `off`.


##Legal stuff

This program is licensed under the GNU GPL version 3.

Copyright 2016 Christoph Hormann
