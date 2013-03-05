#ifndef ZIGNOR_H
#define	ZIGNOR_H

void    RanNormalSetSeedZig(int *piSeed, int cSeed);
double  DRanNormalZig(void);
void    RanNormalSetSeedZigVec(int *piSeed, int cSeed);
double  DRanNormalZigVec(void);
void    RanNormalSetSeedZig32(int *piSeed, int cSeed);
double  DRanNormalZig32(void);
void    RanNormalSetSeedZig32Vec(int *piSeed, int cSeed);
double  DRanNormalZig32Vec(void);

double  DRanQuanNormalZig(void);
double  DRanQuanNormalZigVec(void);
double  DRanQuanNormalZig32(void);
double  DRanQuanNormalZig32Vec(void);

#endif