/* randn_complex_6
  DLW: 27 Feb 2013, 31 Jan 2016
  C-11 complex version

  compile for gcc-5:

gcc -std=gnu11 -Wall -Ofast -msse2 -frename-registers -malign-double -fno-strict-aliasing -DHAVE_SSE2=1 -DDSFMT_MEXP=19937 -o randn_complex_6 dSFMT.c randn_complex_6.c -lm -lrt -lgsl -lgslcblas

*/

#include <complex.h>                // gnu11  complex numbers for C11
#include <math.h>                   // for sine and cosine
#include <time.h>                   // basic C timers, clock, etc    
#include <stdio.h>                  // for keyboard & print IO
#include <stdlib.h>                 // standard C calls
#include <string.h>
#include <sys/time.h>               // for clock_gettime()
#include <assert.h>
#include <gsl/gsl_rng.h>            // needs GNU gsl library under GNU scientific
#include "dSFMT.h"	               // for Mersenne SSE2 generator

dsfmt_t dsfmt;

extern  double  gsl_ran_gaussian_ziggurat (gsl_rng *r, const double sigma);
static  gsl_rng *rng;

/* tabulated values for the height of the Ziggurat levels */
static const double ytab[128] = {
  1.000000000000, 0.963598623011, 0.936280813353, 0.913041104253,
  0.892278506696, 0.873239356919, 0.855496407634, 0.838778928349,
  0.822902083699, 0.807732738234, 0.793171045519, 0.779139726505,
  0.765577436082, 0.752434456248, 0.739669787677, 0.727249120285,
  0.715143377413, 0.703327646455, 0.691780377035, 0.68048276891,
  0.669418297233, 0.65857233912,  0.647931876189, 0.637485254896,
  0.62722199145,  0.617132611532, 0.607208517467, 0.597441877296,
  0.587825531465, 0.578352913803, 0.569017984198, 0.559815170911,
  0.550739320877, 0.541785656682, 0.532949739145, 0.524227434628,
  0.515614886373, 0.507108489253, 0.498704867478, 0.490400854812,
  0.482193476986, 0.47407993601,  0.466057596125, 0.458123971214,
  0.450276713467, 0.442513603171, 0.434832539473, 0.427231532022,
  0.419708693379, 0.41226223212,  0.404890446548, 0.397591718955,
  0.390364510382, 0.383207355816, 0.376118859788, 0.369097692334,
  0.362142585282, 0.355252328834, 0.348425768415, 0.341661801776,
  0.334959376311, 0.328317486588, 0.321735172063, 0.31521151497,
  0.308745638367, 0.302336704338, 0.29598391232,  0.289686497571,
  0.283443729739, 0.27725491156,  0.271119377649, 0.265036493387,
  0.259005653912, 0.253026283183, 0.247097833139, 0.241219782932,
  0.235391638239, 0.229612930649, 0.223883217122, 0.218202079518,
  0.212569124201, 0.206983981709, 0.201446306496, 0.195955776745,
  0.190512094256, 0.185114984406, 0.179764196185, 0.174459502324,
  0.169200699492, 0.1639876086,   0.158820075195, 0.153697969964,
  0.148621189348, 0.143589656295, 0.138603321143, 0.133662162669,
  0.128766189309, 0.123915440582, 0.119109988745, 0.114349940703,
  0.10963544023,  0.104966670533, 0.100343857232, 0.0957672718266,
  0.0912372357329, 0.0867541250127, 0.082318375932, 0.0779304915295,
  0.0735910494266, 0.0693007111742, 0.065060233529, 0.0608704821745,
  0.056732448584,  0.05264727098,   0.0486162607163, 0.0446409359769,
  0.0407230655415, 0.0368647267386, 0.0330683839378, 0.0293369977411,
  0.0256741818288, 0.0220844372634, 0.0185735200577, 0.0151490552854,
  0.0118216532614, 0.00860719483079, 0.00553245272614, 0.00265435214565
};

/* tabulated values for 2^24 times x[i]/x[i+1],
 * used to accept for U*x[i+1]<=x[i] without any floating point operations */
static const unsigned long ktab[128] = {
  0,        12590644, 14272653, 14988939,
  15384584, 15635009, 15807561, 15933577,
  16029594, 16105155, 16166147, 16216399,
  16258508, 16294295, 16325078, 16351831,
  16375291, 16396026, 16414479, 16431002,
  16445880, 16459343, 16471578, 16482744,
  16492970, 16502368, 16511031, 16519039,
  16526459, 16533352, 16539769, 16545755,
  16551348, 16556584, 16561493, 16566101,
  16570433, 16574511, 16578353, 16581977,
  16585398, 16588629, 16591685, 16594575,
  16597311, 16599901, 16602354, 16604679,
  16606881, 16608968, 16610945, 16612818,
  16614592, 16616272, 16617861, 16619363,
  16620782, 16622121, 16623383, 16624570,
  16625685, 16626730, 16627708, 16628619,
  16629465, 16630248, 16630969, 16631628,
  16632228, 16632768, 16633248, 16633671,
  16634034, 16634340, 16634586, 16634774,
  16634903, 16634972, 16634980, 16634926,
  16634810, 16634628, 16634381, 16634066,
  16633680, 16633222, 16632688, 16632075,
  16631380, 16630598, 16629726, 16628757,
  16627686, 16626507, 16625212, 16623794,
  16622243, 16620548, 16618698, 16616679,
  16614476, 16612071, 16609444, 16606571,
  16603425, 16599973, 16596178, 16591995,
  16587369, 16582237, 16576520, 16570120,
  16562917, 16554758, 16545450, 16534739,
  16522287, 16507638, 16490152, 16468907,
  16442518, 16408804, 16364095, 16301683,
  16207738, 16047994, 15704248, 15472926
};

/* tabulated values of 2^{-24}*x[i] */
static const double wtab[128] = {
  1.62318314817e-08, 2.16291505214e-08, 2.54246305087e-08, 2.84579525938e-08,
  3.10340022482e-08, 3.33011726243e-08, 3.53439060345e-08, 3.72152672658e-08,
  3.8950989572e-08, 4.05763964764e-08, 4.21101548915e-08, 4.35664624904e-08,
  4.49563968336e-08, 4.62887864029e-08, 4.75707945735e-08, 4.88083237257e-08,
  5.00063025384e-08, 5.11688950428e-08, 5.22996558616e-08, 5.34016475624e-08,
  5.44775307871e-08, 5.55296344581e-08, 5.65600111659e-08, 5.75704813695e-08,
  5.85626690412e-08, 5.95380306862e-08, 6.04978791776e-08, 6.14434034901e-08,
  6.23756851626e-08, 6.32957121259e-08, 6.42043903937e-08, 6.51025540077e-08,
  6.59909735447e-08, 6.68703634341e-08, 6.77413882848e-08, 6.8604668381e-08,
  6.94607844804e-08, 7.03102820203e-08, 7.11536748229e-08, 7.1991448372e-08,
  7.2824062723e-08, 7.36519550992e-08, 7.44755422158e-08, 7.52952223703e-08,
  7.61113773308e-08, 7.69243740467e-08, 7.77345662086e-08, 7.85422956743e-08,
  7.93478937793e-08, 8.01516825471e-08, 8.09539758128e-08, 8.17550802699e-08,
  8.25552964535e-08, 8.33549196661e-08, 8.41542408569e-08, 8.49535474601e-08,
  8.57531242006e-08, 8.65532538723e-08, 8.73542180955e-08, 8.8156298059e-08,
  8.89597752521e-08, 8.97649321908e-08, 9.05720531451e-08, 9.138142487e-08,
  9.21933373471e-08, 9.30080845407e-08, 9.38259651738e-08, 9.46472835298e-08,
  9.54723502847e-08, 9.63014833769e-08, 9.71350089201e-08, 9.79732621669e-08,
  9.88165885297e-08, 9.96653446693e-08, 1.00519899658e-07, 1.0138063623e-07,
  1.02247952126e-07, 1.03122261554e-07, 1.04003996769e-07, 1.04893609795e-07,
  1.05791574313e-07, 1.06698387725e-07, 1.07614573423e-07, 1.08540683296e-07,
  1.09477300508e-07, 1.1042504257e-07, 1.11384564771e-07, 1.12356564007e-07,
  1.13341783071e-07, 1.14341015475e-07, 1.15355110887e-07, 1.16384981291e-07,
  1.17431607977e-07, 1.18496049514e-07, 1.19579450872e-07, 1.20683053909e-07,
  1.21808209468e-07, 1.2295639141e-07, 1.24129212952e-07, 1.25328445797e-07,
  1.26556042658e-07, 1.27814163916e-07, 1.29105209375e-07, 1.30431856341e-07,
  1.31797105598e-07, 1.3320433736e-07, 1.34657379914e-07, 1.36160594606e-07,
  1.37718982103e-07, 1.39338316679e-07, 1.41025317971e-07, 1.42787873535e-07,
  1.44635331499e-07, 1.4657889173e-07, 1.48632138436e-07, 1.50811780719e-07,
  1.53138707402e-07, 1.55639532047e-07, 1.58348931426e-07, 1.61313325908e-07,
  1.64596952856e-07, 1.68292495203e-07, 1.72541128694e-07, 1.77574279496e-07,
  1.83813550477e-07, 1.92166040885e-07, 2.05295471952e-07, 2.22600839893e-07
};


inline static double gsl_ran_gaussian_zigg (const gsl_rng * r, const double sigma)
{
  double x, y, PARAM_R = 3.44428647676;
  unsigned long int i, j;
  int sign;
  
  const unsigned long int range = r->type->max - r->type->min;
  const unsigned long int offset = r->type->min;

  while (1)
    {
      if (range >= 0xFFFFFFFF)
        {
          unsigned long int k = gsl_rng_get(r) - offset;
          i = (k & 0xFF);
          j = (k >> 8) & 0xFFFFFF;
        }
      else if (range >= 0x00FFFFFF)
        {
          unsigned long int k1 = gsl_rng_get(r) - offset;
          unsigned long int k2 = gsl_rng_get(r) - offset;
          i = (k1 & 0xFF);
          j = (k2 & 0x00FFFFFF);
        }
      else
        {
            i = 256*dsfmt_genrand_uint32(&dsfmt);
            j = 16777216*dsfmt_genrand_uint32(&dsfmt);
        }

      sign = (i & 0x80) ? +1 : -1;
      i &= 0x7f;

      x = j * wtab[i];

      if (j < ktab[i])
        break;

      if (i < 127)
        {
          double y0, y1, U1;
          y0 = ytab[i];
          y1 = ytab[i + 1];
          U1 = genrand_close_open();
          y = y1 + (y0 - y1) * U1;
        }
      else
        {
          double U1, U2;
          U1 = 1.0 - genrand_close_open();
          U2 = genrand_close_open();
          x = PARAM_R - log (U1) / PARAM_R;
          y = exp (-PARAM_R * (x - 0.5 * PARAM_R)) * U2;
        }

      if (y < exp (-0.5 * x * x))
        break;
    }

  return sign * sigma * x;
}

inline static double gauss_ratio (void)
//Joseph Leva algorithm: http://saluc.engr.uconn.edu/refs/crypto/rng/leva92afast.pdf
{
  double u, v, x, y, Q;
  const double s =  0.449871;    /* Constants from Leva */
  const double t = -0.386595;
  const double a =  0.19600;
  const double b =  0.25472;
  const double r1 = 0.27597;
  const double r2 = 0.27846;

  do                            /* This loop is executed 1.369 times on average  */
    {
      u = 1 - dsfmt_genrand_close_open(&dsfmt);
      v = dsfmt_genrand_close_open(&dsfmt) - 0.5;
      v *= 1.7156;
      x = u - s;
      y = fabs (v) - t;
      Q = x * x + y * (a * y - b * x);
    } while (Q >= r1 && (Q > r2 || v * v > -4 * u * u * log (u)));

  return (v / u);       /* Return slope */
}

inline static void randn_polar(long n, complex double r[][n]){
   // Marsaglia polar using both u and v for complex return
   long i, j;
   register double x, y, s
;
   for(i = 0; i < n; i++){
      for(j = 0; j < n; j++){
         for(;;){
            x = 2.*dsfmt_genrand_close_open(&dsfmt) - 1.0;
            y = 2.*dsfmt_genrand_close_open(&dsfmt) - 1.0;
            s = x*x + y*y;
            if (s < 1.0) break;
         }    
         s = sqrt(-2.0*log(s)/s);
         r[i][j] = (x + I*y)*s;
       }
   }
   return;
}

inline static double fast_sin(double(x)){
  double x2; 
  x2 = x*x;                 // defined on -pi to pi by Charles K Garrett, Feb 2012
                            // http://krisgarrett.net/papers/l2approx.pdf

  return ((((((((((+ 2.47852306233493974115e-20*x2  - 8.53932287916564238231e-18)*x2
                   + 2.81875350346861226633e-15)*x2 - 7.64807134493815932275e-13)*x2
                   + 1.60591122567208977895e-10)*x2 - 2.50521116230089813913e-08)*x2
                   + 2.75573193196855760359e-06)*x2 - 1.98412698429672570320e-04)*x2
                   + 8.33333333334987771150e-03)*x2 - 1.66666666666674074058e-01)*x2
                   + 1.00000000000000098216e+00)*x;

/*for 10^-12 error
  return (((((( + 1.35333825545218599272e-10*x2 - 2.47016425480527869032e-08)*x2
                + 2.75322955330449911163e-06)*x2 - 1.98403112669018996690e-04)*x2
                + 8.33331451433080749755e-03)*x2 - 1.66666650437066346286e-01)*x2
                + 9.99999995973569972699e-01)*x;
*/ 
}
 
inline static double fast_cos(double x ){
  double x2;    
  x2 = x*x;                 // defined on -pi to pi by Charles K Garrett, Feb 2012

  return (((((((((+ 3.68396216222400477886e-19*x2 - 1.55289318377801496607e-16)*x2
                  + 4.77840439714556611532e-14)*x2 - 1.14706678499029860238e-11)*x2
                  + 2.08767534780769871595e-09)*x2 - 2.75573191273279748439e-07)*x2
                  + 2.48015873000796780048e-05)*x2 - 1.38888888888779804960e-03)*x2
                  + 4.16666666666665603386e-02)*x2 - 5.00000000000000154115e-01)*x2
                  + 1.00000000000000001607e+00;

/*for 10^-12 error
  return ((((( + 1.73691489450821293670e-09*x2 - 2.71133771940801138503e-07)*x2
               + 2.47734245730930250260e-05)*x2- 1.38879704270452054154e-03)*x2
               + 4.16665243677686230461e-02)*x2 - 4.99999917728614591900e-01)*x2
               + 9.99999992290827491711e-01;
 */ 
}

inline static void randn_box_muller_fast(long n, complex double r[][n])
{
  register double x, y;
  register long i, j;
   
   for(i = 0; i < n; i++){
      for(j = 0; j < n; j++){
         x = 2.*M_PI*dsfmt_genrand_close_open(&dsfmt);
         y = sqrt(-2.*log(dsfmt_genrand_close_open(&dsfmt)));
         r[i][j] = (fast_sin(x) + I*fast_cos(x))*y;
      } 
   } 
  return;
}

inline static void randn_box_muller(long n, complex double r[][n])
{
   long i, j; 
   register double x, y;
  
   for(i = 0; i < n; i++){
      for(j = 0; j < n; j++){  
        x = 2.*M_PI*dsfmt_genrand_close_open(&dsfmt);
        y = sqrt(-2.*log(dsfmt_genrand_close_open(&dsfmt)));
        r[i][j] = (cos(x) + I*sin(x))*y;
      }
    }
  return;
}

inline void SSESqrt_Recip_Times_X( float *pOut, float *pIn )
{
   __m128 in = _mm_load_ss( pIn );
   _mm_store_ss( pOut, _mm_mul_ss( in, _mm_rsqrt_ss( in ) ) );
   // compiles to movss, movaps, rsqrtss, mulss, movss
}

inline static void randn_box_SIMD(long n,complex double r[][n])
{
  
  float yin,yout;
  long i, j;

  register double x;
  for(i = 0; i < n; i++){
      for(j = 0; j < n; j++){
         x = 2.*M_PI*dsfmt_genrand_close_open(&dsfmt);
         yin = (float)(-2*log(dsfmt_genrand_close_open(&dsfmt)));
         SSESqrt_Recip_Times_X(&yout,  &yin );
         r[i][j] = (fast_sin(x)  + I*fast_cos(x))*(double)yout;
      }
  }
  return;
}

double timeDiff( struct timespec *t1, struct timespec *t2){
  double  dt;
    dt = (double)(t2->tv_sec  - t1->tv_sec );
    dt = dt+ (double)(t2->tv_nsec - t1->tv_nsec)/1e+9;
  return dt;
}

int main(void)
{
   struct  timespec t0, t1, t2, t3, t4, t5, t6, t7, t8 ,t9, t10, t11;         // for clock_gettime()
   double u, v;
   long i, j, n, seed;
   const  gsl_rng_type *T;

  n     = 8000;
  seed  = 1234567;

  dsfmt_init_gen_rand(&dsfmt,seed);

/* set up gsl routines      */
  gsl_rng_env_setup();
  T = gsl_rng_default;
  rng = gsl_rng_alloc (T);
  gsl_rng_set(rng, seed);             // seed = 0 default

/* Allocate dynamic, 2-D, complex array with pointers
   from:  http://stackoverflow.com/questions/31511609/multidimensional-arrays-malloc-vs-new   */

   complex double (*r)[n] = calloc(n, sizeof *r);


// Initialize dsfmt with a seed
   dsfmt_init_gen_rand(&dsfmt,seed);

// spin up i7 cpu turbo mode
   randn_box_muller(n, r);
   printf(" u and v = %10.8lf  %10.8lf*I\n", creal(r[0][0]), cimag(r[0][0]));

   clock_gettime(CLOCK_MONOTONIC, &t0);
   for(i = 0; i < n ; i++){
      for(j = 0; j < n ; j++){
         u = gsl_ran_gaussian_zigg (rng, 1.0);
         v = gsl_ran_gaussian_zigg (rng, 1.0);
         r[i][j] = u + I*v;
      }
   }   
   clock_gettime(CLOCK_MONOTONIC, &t1);

   printf(" u and v = %10.8lf  %10.8lf*I\n", creal(r[1][1]), cimag(r[1][1]));

   clock_gettime(CLOCK_MONOTONIC, &t2);
   randn_box_muller(n, r);
   clock_gettime(CLOCK_MONOTONIC, &t3);
   printf(" u and v = %10.8lf  %10.8lf*I\n", creal(r[1][1]), cimag(r[1][1]));

   clock_gettime(CLOCK_MONOTONIC, &t4);
   randn_box_muller_fast(n, r);
   clock_gettime(CLOCK_MONOTONIC, &t5);

   clock_gettime(CLOCK_MONOTONIC, &t6);
   randn_box_SIMD(n, r);
   clock_gettime(CLOCK_MONOTONIC, &t7);
   printf(" u and v = %10.8lf %10.8lf*I\n",creal(r[1][1]), cimag(r[1][1]));
   
   clock_gettime(CLOCK_MONOTONIC, &t8);
   randn_polar(n, r);
   clock_gettime(CLOCK_MONOTONIC, &t9);
   printf(" u and v = %10.8lf  %10.8lf*I\n", creal(r[1][1]), cimag(r[1][1]));

   clock_gettime(CLOCK_MONOTONIC, &t10);
   for(i = 0; i < n; i++){
      for(j = 0; j < n; j++){
         u = gauss_ratio();
         v = gauss_ratio();
         r[i][j] =  u + I*v;
     }
   }
   clock_gettime(CLOCK_MONOTONIC, &t11);
   printf(" u and v = %10.8lf + %10.8lf*I\n", creal(r[1][1]), cimag(r[1][1]));

   printf("\n");

   printf("gsl_zigg_dsfmt                     = %10.8lf (s) \n",  timeDiff(&t0, &t1));
   printf("Box-Muller                         = %10.8lf (s) \n",  timeDiff(&t2, &t3));
   printf("SIMD Box-Muller with fast_sin-cos  = %10.8lf (s) \n",  timeDiff(&t4, &t5));
   printf("Box-Muller manual SIMD             = %10.8lf (s) \n",  timeDiff(&t6, &t7));
   printf("SIMD Marsaglia polar               = %10.8lf (s) \n",  timeDiff(&t8, &t9));
   printf("Time for Leva Gaussian ratio       = %10.8lf (s) \n\n",timeDiff(&t10, &t11));
   
/* Free array memory */
   free(r);
   return 0;
}
