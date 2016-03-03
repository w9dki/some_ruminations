# some_ruminations
some C-11, Fortran , Python and assembly code hacks

These are a comparison of some very fast Gaussian random number generators for complex, double2-D arrays equivalent to the Matlab statement:
  n = 4096;
  E = randn(n, n) + i*randn(n, n) 

or Python3 with numpy / Pylab.  They are meant are expected to be called 10^6 to 10^12 times.

These varients include the Ziggurat, FAST Box-Muller, SIMD Marsaglia Polar, and a Leva Gaussian ratio Gaussian random number algorithms.  All variations use the uniform. Mersenne prime generator dsfmt() and utilize SIMD vectored instructions.

They are all compiled with:   (example for randn_complex_8.c)
 
gcc -std=gnu11 -Wall -Ofast -msse2 -frename-registers -malign-double -fno-strict-aliasing -DHAVE_SSE2=1 -DDSFMT_MEXP=19937 -o randn_complex_8 dSFMT.c randn_complex_8.c -lm -lrt -lgsl -lgslcblas

These are meant to be absolutely as fast as possible.  randn_complex_7.c uses 1_D pointer arithmetic, (~a Fortran equivalence for a 2-D complex, double array) and is not necessarily intended to be particulary easy to understand.  Speed and efficiency are everything.

If you can find a way to make these faster or better.  Let me know.

gcc produces  surprisingly fast machine code in the relatively simple coded 2_D, complex, double versions:  randn_complex_6.c and randn_complex_8.c, whereas randn_complex_5.c and randn_complex_7.c use pointer arithmetic.  In comparison, randn_complex_8.c simply returns a complex pair of pointers ber call and the 2-D complex, double, array allocation and usage occur in main() rather than in the functions themselves.   Pointers to the arrays (calls by location) are used rather than making copies of the huge arrays i.e. calls by value.

The complex gsl ziggurat version with dsfmt() is usually slightly faster, but not dramatically so compared to the simple FAST Box-Muller and Marsaglia Polar algorithms.   gcc  5 gives these timings for the 4090 x 4096 complex double arrays of Gaussian random numbers using 1-D pointer arithmetic on a 2-D complex double array  (normally, the fastest version for a 4096 x 4096 double complex array in gcc, C11

$ ./randn_complex_7
 u and v = -0.27161188  0.12163298*I
 u and v = -1.50471852  -0.64422578*I
 u and v = -0.41339744  -0.82415066*I
 u and v = 0.28028943  -0.25109037*I
 u and v = 1.12512160  -0.56908906*I
 u and v = -1.18069878  0.61920043*I
 
gsl_zigg_dsfmt                     = 0.43400046 (s)
Box-Muller                         = 1.37106105 (s) 
SIMD Box-Muller with fast_sin-cos  = 0.63448656 (s)
Box-Muller manual SIMD             = 0.65211192 (s)
SIMD Marsaglia polar               = 0.62771889 (s)
Time for Leva Gaussian ratio       = 0.40401675 (s)

**************

Dividing by 4096 x 4096 x 2 gives 12.8 ns per Gaussian random number for the gsl, Ziggurat code using dsfmt() and 18 ns per Gaussian random number with the SIMD Marsaglia polar varient.


DLW/W9DKI
