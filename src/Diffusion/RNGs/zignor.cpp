/*==========================================================================
 *  This code is Copyright (C) 2005, Jurgen A. Doornik.
 *  Permission to use this code for non-commercial purposes
 *  is hereby given, provided proper reference is made to:
 *		Doornik, J.A. (2005), "An Improved Ziggurat Method to Generate Normal
 *          Random Samples", mimeo, Nuffield College, University of Oxford,
 *			and www.doornik.com/research.
 *		or the published version when available.
 *	This reference is still required when using modified versions of the code.
 *  This notice should be maintained in modified versions of the code.
 *	No warranty is given regarding the correctness of this code.
 *==========================================================================*/
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "zigrandom.h"
#include "zignor.h"


/*------------------------------ General Ziggurat --------------------------*/
static double DRanNormalTail(double dMin, int iNegative)
{
	double x, y;
	do
	{	x = log(DRanU()) / dMin;
		y = log(DRanU());
	} while (-2 * y < x * x);
	return iNegative ? x - dMin : dMin - x;
}

#define ZIGNOR_C 128			       /* number of blocks */
#define ZIGNOR_R 3.442619855899	/* start of the right tail */
				   /* (R * phi(R) + Pr(X>=R)) * sqrt(2\pi) */
#define ZIGNOR_V 9.91256303526217e-3

/* s_adZigX holds coordinates, such that each rectangle has*/
/* same area; s_adZigR holds s_adZigX[i + 1] / s_adZigX[i] */
static double s_adZigX[ZIGNOR_C + 1], s_adZigR[ZIGNOR_C];

static void zigNorInit(int iC, double dR, double dV)
{
	int i;	double f;
	
	f = exp(-0.5 * dR * dR);
	s_adZigX[0] = dV / f; /* [0] is bottom block: V / f(R) */
	s_adZigX[1] = dR;
	s_adZigX[iC] = 0;

	for (i = 2; i < iC; ++i)
	{
		s_adZigX[i] = sqrt(-2 * log(dV / s_adZigX[i - 1] + f));
		f = exp(-0.5 * s_adZigX[i] * s_adZigX[i]);
	}
	for (i = 0; i < iC; ++i)
		s_adZigR[i] = s_adZigX[i + 1] / s_adZigX[i];
}
double  DRanNormalZig(void)
{
	unsigned int i;
	double x, u, f0, f1;
	
	for (;;)
	{
		u = 2 * DRanU() - 1;
		i = IRanU() & 0x7F;
		/* first try the rectangular boxes */
		if (fabs(u) < s_adZigR[i])		 
			return u * s_adZigX[i];
		/* bottom box: sample from the tail */
		if (i == 0)						
			return DRanNormalTail(ZIGNOR_R, u < 0);
		/* is this a sample from the wedges? */
		x = u * s_adZigX[i];		   
		f0 = exp(-0.5 * (s_adZigX[i] * s_adZigX[i] - x * x) );
		f1 = exp(-0.5 * (s_adZigX[i+1] * s_adZigX[i+1] - x * x) );
      	if (f1 + DRanU() * (f0 - f1) < 1.0)
			return x;
	}
}

#define ZIGNOR_STORE 64 * 4
static unsigned int s_auiZigTmp[ZIGNOR_STORE / 4];
static unsigned int s_auiZigBox[ZIGNOR_STORE];
static double s_adZigRan[ZIGNOR_STORE + ZIGNOR_STORE / 4];
static int s_cZigStored = 0;

double  DRanNormalZigVec(void)
{
	unsigned int i, j, k;
	double x, u, f0, f1;
	
	for (;;)
	{
		if (s_cZigStored == 0)
		{
			RanVecIntU(s_auiZigTmp, ZIGNOR_STORE / 4);
			RanVecU(s_adZigRan, ZIGNOR_STORE);
			for (j = k = 0; j < ZIGNOR_STORE; j += 4, ++k)
			{
				i = s_auiZigTmp[k];	s_auiZigBox[j + 0] = i & 0x7F;
				i >>= 8;			s_auiZigBox[j + 1] = i & 0x7F;
				i >>= 8;			s_auiZigBox[j + 2] = i & 0x7F;
				i >>= 8;			s_auiZigBox[j + 3] = i & 0x7F;
				s_adZigRan[j + 0] = 2 * s_adZigRan[j + 0] - 1;
				s_adZigRan[j + 1] = 2 * s_adZigRan[j + 1] - 1;
				s_adZigRan[j + 2] = 2 * s_adZigRan[j + 2] - 1;
				s_adZigRan[j + 3] = 2 * s_adZigRan[j + 3] - 1;
			}
			s_cZigStored = j;
		}
		--s_cZigStored;

		u = s_adZigRan[s_cZigStored];
		i = s_auiZigBox[s_cZigStored];
		
		if (fabs(u) < s_adZigR[i])		 /* first try the rectangular boxes */
			return u * s_adZigX[i];

		if (i == 0)						/* bottom box: sample from the tail */
			return DRanNormalTail(ZIGNOR_R, u < 0);

		x = u * s_adZigX[i];		   /* is this a sample from the wedges? */
		f0 = exp(-0.5 * (s_adZigX[i] * s_adZigX[i] - x * x) );
		f1 = exp(-0.5 * (s_adZigX[i + 1] * s_adZigX[i + 1] - x * x) );
      	if (f1 + DRanU() * (f0 - f1) < 1.0)
			return x;
	}
}

void  RanNormalSetSeedZig(int *piSeed, int cSeed)
{
	zigNorInit(ZIGNOR_C, ZIGNOR_R, ZIGNOR_V);
	RanSetSeed(piSeed, cSeed);
}
void  RanNormalSetSeedZigVec(int *piSeed, int cSeed)
{
	s_cZigStored = 0;
	RanNormalSetSeedZig(piSeed, cSeed);
}
/*--------------------------- END General Ziggurat -------------------------*/

/*------------------------------ Integer Ziggurat --------------------------*/
#define ZIGNOR_INVM	M_RAN_INVM32

static unsigned int s_aiZigRm[ZIGNOR_C];
static double s_adZigXm[ZIGNOR_C + 1];

static void zig32NorInit(int iC, double dR, double dV)
{
	int i;  double f, m31 = ZIGNOR_INVM * 2;

	f = exp(-0.5 * dR * dR);
	s_adZigXm[0] = dV / f; /* [0] is bottom block: V / f(R) */
	s_adZigXm[1] = dR;
	s_adZigXm[iC] = 0;

	for (i = 2; i < iC; ++i)
	{
		s_adZigXm[i] = sqrt(-2 * log(dV / s_adZigXm[i - 1] + f));
		f = exp(-0.5 * s_adZigXm[i] * s_adZigXm[i]);
	}
	/* compute ratio and implement scaling */
	for (i = 0; i < iC; ++i)
		s_aiZigRm[i] = (unsigned int)
			( (s_adZigXm[i + 1] / s_adZigXm[i]) / m31 );
	for (i = 0; i <= iC; ++i)
		s_adZigXm[i] *= m31;
}
double  DRanNormalZig32(void)
{
	unsigned int i;
	int u;
	double x, y, f0, f1;
	
	for (;;)
	{
		u = (int)IRanU();
		i = IRanU() & 0x7F;
		if ((unsigned int)abs(u) < s_aiZigRm[i])/* first try the rectangles */
			return u * s_adZigXm[i];
		
		if (i == 0)									/* sample from the tail */
			return DRanNormalTail(ZIGNOR_R, u < 0);

		x = u * s_adZigXm[i];		   /* is this a sample from the wedges? */
		y = 0.5 * s_adZigXm[i] / ZIGNOR_INVM;      f0 = exp(-0.5 * (y * y - x * x) );
		y = 0.5 * s_adZigXm[i + 1] / ZIGNOR_INVM;  f1 = exp(-0.5 * (y * y - x * x) );
      	if (f1 + IRanU() * ZIGNOR_INVM * (f0 - f1) < 1.0)
			return x;
	}
}
#define ZIGNOR32_STORE 64 * 4
static unsigned int s_auiZig32Ran[ZIGNOR32_STORE];
static unsigned int s_auiZig32Box[ZIGNOR32_STORE];
static int s_cZig32Stored = 0;

double  DRanNormalZig32Vec(void)
{
	unsigned int i, j, k;
	int u;
	double x, y, f0, f1;
	
	for (;;)
	{
		if (s_cZig32Stored == 0)
		{
			RanVecIntU(s_auiZig32Ran, ZIGNOR32_STORE / 4);
			for (j = k = 0; j < ZIGNOR32_STORE; j += 4, ++k)
			{
				i = s_auiZig32Ran[k];	s_auiZig32Box[j+0] = i & 0x7F;
				i >>= 8;				s_auiZig32Box[j+1] = i & 0x7F;
				i >>= 8;				s_auiZig32Box[j+2] = i & 0x7F;
				i >>= 8;				s_auiZig32Box[j+3] = i & 0x7F;
			}
			RanVecIntU(s_auiZig32Ran, ZIGNOR32_STORE);
			s_cZig32Stored = j;
		}
		--s_cZig32Stored;

		u = (int)s_auiZig32Ran[s_cZig32Stored];
		i = s_auiZig32Box[s_cZig32Stored];
		/* first try the rectangles */
		if ((unsigned int)abs(u) < s_aiZigRm[i])
			return u * s_adZigXm[i];
		/* bottom box: sample from the tail */
		if (i == 0)									
			return DRanNormalTail(ZIGNOR_R, u < 0);
		/* is this a sample from the wedges? */
		x = u * s_adZigXm[i];
		y = 0.5 * s_adZigXm[i] / ZIGNOR_INVM;
		f0 = exp(-0.5 * (y * y - x * x) );
		y = 0.5 * s_adZigXm[i + 1] / ZIGNOR_INVM;
		f1 = exp(-0.5 * (y * y - x * x) );
      	if (f1 + IRanU() * ZIGNOR_INVM * (f0 - f1) < 1.0)
			return x;
	}
}
void  RanNormalSetSeedZig32(int *piSeed, int cSeed)
{
	zig32NorInit(ZIGNOR_C, ZIGNOR_R, ZIGNOR_V);
	RanSetSeed(piSeed, cSeed);
}
void  RanNormalSetSeedZig32Vec(int *piSeed, int cSeed)
{
	s_cZig32Stored = 0;
	RanNormalSetSeedZig32(piSeed, cSeed);
}
/*--------------------------- END Integer Ziggurat -------------------------*/

/*--------------------------- functions for testing ------------------------*/
double  DRanQuanNormalZig(void)
{
	return DProbNormal(DRanNormalZig());
}
double  DRanQuanNormalZigVec(void)
{
	return DProbNormal(DRanNormalZigVec());
}
double  DRanQuanNormalZig32(void)
{
	return DProbNormal(DRanNormalZig32());
}
double  DRanQuanNormalZig32Vec(void)
{
	return DProbNormal(DRanNormalZig32Vec());
}
/*------------------------- END functions for testing ----------------------*/
