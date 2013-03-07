/* 
 * File:   QMCheaders.h
 * Author: jorgehog
 *
 * Created on 26. juni 2012, 15:19
 */

#ifndef QMCHEADERS_H
#define	QMCHEADERS_H

//#define ARMA_NO_DEBUG
//#define OMP_ON
#define MPI_ON

#ifdef OMP_ON
#include <omp.h>
#endif
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



#include <armadillo>
#include <boost/lexical_cast.hpp>
#include <stdlib.h>
#include <iostream> 
#include <fstream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <sys/time.h>

#define TOSTR boost::lexical_cast<std::string>


struct DMCparams {
    int n_c, therm, n_b, n_w;
    double dt;
};

struct VMCparams {
    int n_c;
    double dt;
};

struct VariationalParams {
    double alpha, beta;
};

struct ParParams {
    bool parallel;
    bool is_master;
    int n_nodes;
    int node;
};

struct GeneralParams {
    int n_p, dim;
    seed_type random_seed;

    double systemConstant;

    bool doMIN;
    bool doVMC;
    bool doDMC;

    bool do_blocking;

    bool use_jastrow;
    bool use_coulomb;

    std::string system;
    std::string sampling;
    
    std::string runpath;

};

struct MinimizerParams {
    double max_step;
    double f_max;
    double f_min;
    double omega;
    double A;
    double a;
    int SGDsamples;
    int n_w;
    int therm;
    int n_c;
    int n_c_SGD;
    arma::rowvec alpha;
    arma::rowvec beta;

};

struct OutputParams {
    bool dist_out;
    bool dmc_out;
    bool ASGD_out;
};

class Orbitals;
class Potential;
class System;
class Sampling;
class Jastrow;

struct SystemObjects {
    Orbitals* SP_basis;
    Potential* onebody_pot;
    System* SYSTEM;
    Sampling* sample_method;
    Jastrow* jastrow;

};

class STDOUT {
public:
    
    virtual void cout(std::stringstream & a) {
        std::cout << a.str() << std::endl;
        a.str(std::string());
        a.clear();
    }
};

class NO_STDOUT : public STDOUT {
public:

    virtual void cout(std::stringstream & a) {
        a.str(std::string());
        a.clear();
    }
};

#include "Walker/Walker.h"

class Minimizer;

class OutputHandler;

class ErrorEstimator;

class QMC;

#include "BasisFunctions/BasisFunctions.h"
#include "BasisFunctions/HarmonicOscillator/HarmonicOscillator.h"
#include "BasisFunctions/hydrogenic/hydrogenic.h"

#include "Orbitals/Orbitals.h"
#include "Orbitals/AlphaHarmonicOscillator/AlphaHarmonicOscillator.h"
#include "Orbitals/hydrogenicOrbitals/hydrogenicOrbitals.h"
#include "Orbitals/ExpandedBasis/ExpandedBasis.h"

#include "Potential/Potential.h"
#include "Potential/Coulomb/Coulomb.h"
#include "Potential/Harmonic_osc/Harmonic_osc.h"
#include "Potential/AtomCore/AtomCore.h"

#include "Jastrow/Jastrow.h"
#include "Jastrow/Pade_Jastrow/Pade_Jastrow.h"
#include "Jastrow/No_Jastrow/No_Jastrow.h"

#include "System/System.h"
#include "System/Fermions/Fermions.h"

#include "QMC/QMC.h"
#include "QMC/VMC/VMC.h"
#include "QMC/DMC/DMC.h"

#include "Diffusion/Diffusion.h"
#include "Diffusion/Simple/Simple.h"
#include "Diffusion/Fokker_Planck/Fokker_Planck.h"

#include "Sampling/Sampling.h"
#include "Sampling/Brute_Force/Brute_Force.h"
#include "Sampling/Importance/Importance.h"

#include "ErrorEstimator/ErrorEstimator.h"
#include "ErrorEstimator/Blocking/Blocking.h"
#include "ErrorEstimator/SimpleVar/SimpleVar.h"

#include "Minimizer/Minimizer.h"
#include "Minimizer/ASGD/ASGD.h"

#include "OutputHandler/OutputHandler.h"
#include "OutputHandler/Distribution/Distribution.h"
#include "OutputHandler/stdoutDMC/stdoutDMC.h"
#include "OutputHandler/stdoutASGD/stdoutASGD.h"

#endif	/* QMCHEADERS_H */

