/* 
 * File:   Diffusion.cpp
 * Author: jorgehog
 * 
 * Created on 16. april 2012, 14:03
 */

#include "../QMCheaders.h"


Diffusion::Diffusion(int n_p, int dim, double timestep, long random_seed, double D) {
    this->n_p = n_p;
    this->dim = dim;
    this->timestep = timestep;
    this->random_seed = random_seed;
    this->D = D;
    this->std = sqrt(2 * D * timestep);
}

double Diffusion::gaussian_deviate(long * idum) {
    static int iset = 0;
    static double gset;
    double fac, rsq, v1, v2;

    if (idum < 0) iset = 0;
    if (iset == 0) {
        do {
            v1 = 2. * ran2(idum) - 1.0;
            v2 = 2. * ran2(idum) - 1.0;
            rsq = v1 * v1 + v2*v2;
        } while (rsq >= 1.0 || rsq == 0.);
        fac = sqrt(-2. * log(rsq) / rsq);
        gset = v1*fac;
        iset = 1;
        return v2*fac;
    } else {
        iset = 0;
        return gset;
    }
}


#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double Diffusion::ran2(long *idum)
{
   int            j;
   long           k;
   static long    idum2 = 123456789;
   static long    iy=0;
   static long    iv[NTAB];
   double         temp;

   if(*idum <= 0) {
      if(-(*idum) < 1) *idum = 1;
      else             *idum = -(*idum);
      idum2 = (*idum);
      for(j = NTAB + 7; j >= 0; j--) {
         k     = (*idum)/IQ1;
	 *idum = IA1*(*idum - k*IQ1) - k*IR1;
	 if(*idum < 0) *idum +=  IM1;
	 if(j < NTAB)  iv[j]  = *idum;
      }
      iy=iv[0];
   }
   k     = (*idum)/IQ1;
   *idum = IA1*(*idum - k*IQ1) - k*IR1;
   if(*idum < 0) *idum += IM1;
   k     = idum2/IQ2;
   idum2 = IA2*(idum2 - k*IQ2) - k*IR2;
   if(idum2 < 0) idum2 += IM2;
   j     = iy/NDIV;
   iy    = iv[j] - idum2;
   iv[j] = *idum;
   if(iy < 1) iy += IMM1;
   if((temp = AM*iy) > RNMX) return RNMX;
   else return temp;
}

#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX