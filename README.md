# some_ruminations
some C-11, Fortran , Python and assembly code hacks

These first are very fast random number generators for complex, double random number generators equivalent to
randn() + i*randn()    in Matlab or Python3 with numpy / Pylab.  They are meant are meant to be called 10^6 to 10^12 times.

These varients include the Ziggurat, FAST Box-Muller, SIMD Marsaglia Polar, and a Leva Gaussian ratio algorithm.  All variatioins use dsfmt() and SIMD vectored instructions.  They are compiled with:   (example for randn_complex_8.c)
 
gcc -std=gnu11 -Wall -Ofast -msse2 -frename-registers -malign-double -fno-strict-aliasing -DHAVE_SSE2=1 -DDSFMT_MEXP=19937 -o randn_complex_8 dSFMT.c randn_complex_8.c -lm -lrt -lgsl -lgslcblas

These are meant to be absolutely as fast as possible.  randn_compolex_7.c uses 1_D pointer arithmetic, (~a Fortran equivalence for a 2-D complex, double array) and are not intended to be particulary easy to understand.  Speed and efficiency are everything.

If you can find a way to make these faster or better.  Let me know.
\
gcc produces  surprisingly fast mavhine code in the relatively simple 2_D, complex, double versions  randn_complex_6.c and randn_complex_8.c versions, _5.c and _7.c use pointer arithmetic, and _8.c simply returns a complex pairof pointers  rather than a full, 2_D complex array with pointers.

The complex gsl ziggurat version with dsfmt() is usually slightly faster, but not dramatically so compared to the simple FAST Box-Muller and Marsaglia Polar algorithms

DLW/W9DKI
