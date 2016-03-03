# some_ruminations
some C-11, Fortran , Python and assembly code hacks

These first are very fast random number generators for complex, double random number generators equivalent to
randn() + i*randn()    in Matlab or Python3 with numpy / matlib.  They are meant are meant to be called 10^6 to 10^12 times.

These varients include the Ziggurat, FAST Box-Muller, SIMD Marsaglia Polar, and a Leva Gaussian ratio algorithm.  All variatioins use dsfmt() and SIMD vectored instructions.  They are compiled with:   (example for randn_complex_8.c)
 
gcc -std=gnu11 -Wall -Ofast -msse2 -frename-registers -malign-double -fno-strict-aliasing -DHAVE_SSE2=1 -DDSFMT_MEXP=19937 -o randn_complex_8 dSFMT.c randn_complex_8.c -lm -lrt -lgsl -lgslcblas

These are meant to be absolutely as fast as possible.  randn_complex_7.c uses a Fortran like equivalence ( i.e. somewhat like a C-union) accessing the 2-D, complex, double array as a 1-D real array and pointer arithmetic.  Note, these are not intended to be particulary easy to understand.  Speed and efficiency are absolutely everything.  I am counting machine instructions with a gcc -S assembly code dump.

If you can find a way to make these faster or better in any way.  Let me know.   Surprisingly to me, gcc produces  very fast machine code in the relatively simple 2-D, complex, double, versions  randn_complex_6.c and randn_complex_8.c,  My hand coded _5.c and _7.c versuse use pointer arithmetic,  Even re surprisingly, randn_complex _8.c simply returns a only a pair of pointers  for the complex variable, rather than a full, 2_D complex array with pointers and is almost as fast as the hand coded, pointer arithmetic, register variable versions.

The complex gsl ziggurat version with dsfmt() is usually the fastest, but not by much compared to the simple fast -Box-Muller and Marsaglia Polar algorithms.   Note the use of the fast_cosine() and fast_sin() algorithms of Charles K Garrett for a factor of ~3 speedup of the Box-Muller algorithm.

DLW/W9DKI
Mar 2, 2016
