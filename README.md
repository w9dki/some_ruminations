# some_ruminations
some C-11, Fortran , Python and assembly code hacks

These first are very fast random number generators for complex, double random number generators equivalent to
randn() + i*randn()

in Matlab or Python3 with numpy / matlib.  They are meant are meant to be called 10^6 to 10^12 times.

These varients include the Ziggurat, Box Muller,SIMD Marsaglia polar, and Leva Gaussian ratio and all use dsfmt() and SIMD vectored instructions.  They are compiled with:   (example for randn_complex_8.c)
 
gcc -std=gnu11 -Wall -Ofast -msse2 -frename-registers -malign-double -fno-strict-aliasing -DHAVE_SSE2=1 -DDSFMT_MEXP=19937 -o randn_complex_8 dSFMT.c randn_complex_8.c -lm -lrt -lgsl -lgslcblas

These are meant to be absolutely as fast as possible.  randn_compolex_7.c uses 1_D pointer arithmetic, (~a Fortran equivalence for a 2-D complex, double array) and are not intended to be particulary easy to understand.  Speed and efficiency are everything.

If you can find a way to make these faster or better.  Let me know.

DLW/W9DKI
