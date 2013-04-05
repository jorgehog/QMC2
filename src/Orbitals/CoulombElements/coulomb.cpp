#include <math.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "coulomb.h"


// compute one non-antisymmetric part of the Coulomb matrix element: <12|V|34> [see Asinimovas and Matulis (1997) for analytical expression]
// does NOT compute <12|V|34>as.
// Be careful, this function doesn't account for spin. 
//
// ex.: exchangeT= \delta(ms1,ms4) * \delta(ms2,ms3) * anisimovas(n1,ml1,n2,ml2,n3,ml3,n4,ml4);
//      directT  = \delta(ms1,ms3) * \delta(ms2,ms4) * anisimovas(n1,ml1,n2,ml2,n4,ml4,n3,ml3);
//
double coulomb (const int n1, const int m1, const int n2, const int m2, const int n3, const int m3, const int n4, const int m4)
{
  double coulombint=0.0;
  if(m1+m2 == m3+m4)
    {
      double temp;
      int lambda;

      int gamma1=0;
      int gamma2=0;
      int gamma3=0;
      int gamma4=0;
      int G=0;
      for(int j1=0;j1<=n1;j1++)
	for(int j2=0;j2<=n2;j2++)
	  for(int j3=0;j3<=n3;j3++)
	    for(int j4=0;j4<=n4;j4++)
	      {
		gamma1 = (int) (j1 +j4 + 0.5*(abs(m1)+m1) + 0.5*(abs(m4)-m4));
		gamma2 = (int) (j2 +j3 + 0.5*(abs(m2)+m2) + 0.5*(abs(m3)-m3));
		gamma3 = (int) (j2 +j3 + 0.5*(abs(m3)+m3) + 0.5*(abs(m2)-m2));
		gamma4 = (int) (j1 +j4 + 0.5*(abs(m4)+m4) + 0.5*(abs(m1)-m1));

		G=gamma1+gamma2+gamma3+gamma4;

		double Lgratio1 = LogRatio1(j1,j2,j3,j4);
		double Lgproduct2 = LogProduct2(n1,m1,n2,m2,n3,m3,n4,m4,j1,j2,j3,j4);
		double Lgratio2 = LogRatio2(G);

		temp=0.0;
		lambda=0;
		for(int l1=0;l1<=gamma1;l1++)
		  for(int l2=0;l2<=gamma2;l2++)
		    for(int l3=0;l3<=gamma3;l3++)
		      for(int l4=0;l4<=gamma4;l4++)
			{
			  lambda=l1+l2+l3+l4;
			  if( (l1+l2)==(l3+l4) )
			    temp += minusPower(gamma2+gamma3-l2-l3)*exp(LogProduct3(l1,l2,l3,l4,gamma1,gamma2,gamma3,gamma4)+lgamma(1+lambda*0.5)+lgamma((G-lambda+1)*0.5));
			}


		coulombint += minusPower(j1+j2+j3+j4) *exp(Lgratio1 + Lgproduct2 + Lgratio2) *temp;

	      }
      coulombint *= Product1(n1,m1,n2,m2,n3,m3,n4,m4);
    }
  return coulombint;
}


// Computes the first ratio in the Asinimovas expression
double LogRatio1(const int j1,const int j2,const int j3,const int j4)
{
  double temp = -LogFac(j1)-LogFac(j2)-LogFac(j3)-LogFac(j4);
  return temp;
} 

// Computes the 2nd ratio in the Asinimovas expression
double LogRatio2(const int G)
{
  double temp = -1*(G+1)*0.5*log(2);
  return temp;
} 


// Computes the log of the 2nd product in the Asinimovas expression
double LogProduct2(const int n1,const int m1,const int n2,const int m2,const int n3,const int m3,const int n4,const int m4,const int j1,const int j2,const int j3,const int j4)
{
  double temp = LogFac(n1+abs(m1))+LogFac(n2+abs(m2))+LogFac(n3+abs(m3))+LogFac(n4+abs(m4))-LogFac(n1-j1)-LogFac(n2-j2)-LogFac(n3-j3)-LogFac(n4-j4)-LogFac(j1+abs(m1))-LogFac(j2+abs(m2))-LogFac(j3+abs(m3))-LogFac(j4+abs(m4));
  return temp;
}


// Computes the log of the 3rd product in the Asinimovas expression
double LogProduct3(const int l1,const int l2,const int l3,const int l4,const int gamma1,const int gamma2,const int gamma3,const int gamma4)
{
  double temp = LogFac(gamma1)+LogFac(gamma2)+LogFac(gamma3)+LogFac(gamma4)-LogFac(l1)-LogFac(l2)-LogFac(l3)-LogFac(l4)-LogFac(gamma1-l1)-LogFac(gamma2-l2)-LogFac(gamma3-l3)-LogFac(gamma4-l4);
  return temp;
}


// computes (-1)^k
int minusPower(const int k)
{
  int temp=abs(k%2)*-2+1;// gives 1 if k is even, gives -1 if k is odd
  return temp;
}


// computes log(n!)
double LogFac (const int n)
{
  if(n>170)
    {
      printf("#### Too big integer in LogFac(n)!!!!#### \n\n");
      exit(1);
    }
  double temp=0;
  for(int i=2; i<n+1;i++)
    temp+=log(i);
  return temp;
}

// computes first product of indices in the Anisimovas/Matulis expression
double Product1 (const int n1, const int ml1, const int n2, const int ml2, const int n3, const int ml3, const int n4, const int ml4)
{
  double temp=0;
  temp = LogFac(n1)+LogFac(n2)+LogFac(n3)+LogFac(n4)-LogFac(n1+abs(ml1))-LogFac(n2+abs(ml2))-LogFac(n3+abs(ml3))-LogFac(n4+abs(ml4));
  temp*=0.5;
  return exp(temp);
}

// lgamma.cpp -- log gamma function of real argument.
//      Algorithms and coefficient values from "Computation of Special
//      Functions", Zhang and Jin, John Wiley and Sons, 1996.
//
//  (C) 2003, C. Bond. All rights reserved.
//
//  Returns log(gamma) of real argument.
//  NOTE: Returns 1e308 if argument is 0 or negative.
//  taken on the web [http://www.crbond.com/math.htm]
double lgamma(double x)
{
#if CHECK
  if(x>=171)
    {
      printf("#### Too big integer in lgamma(x) to give accurate result!!!!#### \n\n");
      exit(1);
    }
#endif  
  
  double x0,x2,xp,gl,gl0;
  int n,k;
  static double a[] = {
    8.333333333333333e-02,
    -2.777777777777778e-03,
    7.936507936507937e-04,
    -5.952380952380952e-04,
    8.417508417508418e-04,
    -1.917526917526918e-03,
    6.410256410256410e-03,
    -2.955065359477124e-02,
    1.796443723688307e-01,
    -1.39243221690590};
    
  x0 = x;
  if (x <= 0.0) return 1e308;
  else if ((x == 1.0) || (x == 2.0)) return 0.0;
  else if (x <= 7.0) {
    n = (int)(7-x);
    x0 = x+n;
  }
  x2 = 1.0/(x0*x0);
  xp = 2.0*M_PI;
  gl0 = a[9];
  for (k=8;k>=0;k--) {
    gl0 = gl0*x2 + a[k];
  }
  gl = gl0/x0+0.5*log(xp)+(x0-0.5)*log(x0)-x0;
  if (x <= 7.0) {
    for (k=1;k<=n;k++) {
      gl -= log(x0-1.0);
      x0 -= 1.0;
    }
  }
  return gl;
}// end of lgamma function
