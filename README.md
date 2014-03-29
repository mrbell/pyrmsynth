`pyrmsynth` - Python based RM Synthesis code including RMCLEAN
==============================================================

*Current version:* 1.2.1
*Updated on:* 2014-03-29

`pyrmsynth` performs RM synthesis, either simply by Fourier transformation 
(to produce a dirty image) or using the RMCLEAN method as described by 
Heald, et al. (2009).  It uses FFTs for the Fourier inversion and, as far as I 
know, this is the only RM synthesis software around that does this. The Numpy 
FFTs are themselves quite fast, but in order to use them, the data first need to 
be placed on a regularly spaced lambda^2 grid. For this, the data are "gridded" 
by convolution with a Kaiser-Bessel Window function and sampling at regular 
intervals, as described in e.g. Beatty, Nishimura, and Pauly 
(IEEE Transactions in Med. Imaging, Vol 24, No. 6, 2005).

The gridding procedure, which requires a convolution, is quite slow when 
implemented in pure Python, so I have re-implemented the gridding routines using
Cython, which converts python-eque code into C code that can be compiled and 
imported into Python.

The result is a package that performs fast RM Synthesis and RM CLEAN imaging
while still providing the flexibility of a Python interface.

pyrmsynth contains an application `rmsynthesis.py` for processing the lines of 
sight in a "stack" of sky images, i.e. many polarized sky images at different 
frequencies. This application was written with [LOFAR](http://www.lofar.org) 
processing in mind, but should be generally useful as long as your images are
provided as a set of FITS files generaged by CASA (or something compatible).

pyrmsynth also can be used as a more generic library for writing your own 
RM synthesis applications. The `rm_tools` sub-package contains efficient classes
for RM synthesis and RM CLEAN. You can easily write your own scripts to do file
I/O and use this package to do the actual RM synthesis computations in an
efficient manner.

For more information, please refer to the [pyrmsynth Wiki](https://github.com/mrbell/gfft/wiki).

`pyrmsynth` is licensed under the [GPLv3](http://www.gnu.org/licenses/gpl.html).
