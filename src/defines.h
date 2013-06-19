#ifndef CFG_H
#define CFG_H 

#define ARMA_NO_DEBUG
#define MPI_ON

#ifdef MPI_ON
#include <mpi.h>
#endif

//#define RNG_NUMREC
#define RNG_ZIG

#ifdef RNG_ZIG
#ifdef RNG_NUMREC
#undef RNG_NUMREC
#endif
typedef int seed_type;
#endif

#ifdef RNG_NUMREC
typedef long seed_type;
#endif

#ifndef TOSTR
#include <boost/lexical_cast.hpp>
#define TOSTR boost::lexical_cast<std::string>
#endif

#endif	/* CFG_H */