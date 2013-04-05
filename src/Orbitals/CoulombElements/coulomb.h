// compute a part of the Coulomb matrix element: <12|V|34>
double coulomb(const int, const int, const int, const int, const int, const int, const int, const int);

//compute (-1)^k
int minusPower(const int);

// computes log(n!)
double LogFac (const int);

// Computes the first ratio in the Asinimovas expression
double LogRatio1(const int,const int,const int,const int);

// Computes the 2nd ratio in the Asinimovas expression
double LogRatio2(const int);

// computes first product of indices in the Anisimovas/Matulis expression
double Product1 (const int, const int, const int, const int, const int, const int, const int, const int);

// Computes the log of the 2nd product in the Asinimovas expression
double LogProduct2(const int,const int,const int,const int,const int,const int,const int,const int,const int,const int,const int,const int);

// Computes the log of the 3rd product in the Asinimovas expression
double LogProduct3(const int,const int,const int,const int,const int,const int,const int,const int);

// The function lgamma() computes the logarithm of the gamma function of real argument x
double lgamma(double x);
