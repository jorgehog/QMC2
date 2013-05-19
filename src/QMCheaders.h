/* 
 * File:   QMCheaders.h
 * Author: jorgehog
 *
 * Created on 26. juni 2012, 15:19
 */

#ifndef QMCHEADERS_H
#define	QMCHEADERS_H

//#define ARMA_NO_DEBUG
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


//! Struct used to initialize DMC parameters.
struct DMCparams {
    int n_c; //!< The number of cycles
    int therm; //!< Thermalization cycles
    int n_b; //!< Number of block samples pr. walker pr. cycle
    int n_w; //!< Number of walkers.
    double dt; //!< Time step
};

//! Struct used to initialize VMC parameters.
struct VMCparams {
    int n_c; //!< The number of cycles.
    double dt; //!< The time step.
};

//! Struct used to initialize the varational parameters.
struct VariationalParams {
    double alpha; //!< The spatial variational parameter.
    double beta; //!< The Jastrow variational parameter.
};

//! Struct used to initialize parallelization parameters.
struct ParParams {
    bool is_master; //!< True for the master node.
    bool parallel; //!< True if n_nodes > 1.
    int node; //!< The process' rank.
    int n_nodes; //!< The total number of processes.
};

//! Struct used to initialize general parameters.
struct GeneralParams {
    int n_p; //!< The number of particles.
    int dim; //!< The dimension.
    seed_type random_seed; //!< The random number generator's seed.
    
    //! The constant used in systems
    /*!
     * e.g. charge for atoms and oscillator frequency for quantum dots.
     */
    double systemConstant; 

    double R; //! Center of mass coordinate for diatmic systems.
    
    bool doMIN; 
    bool doVMC;
    bool doDMC;

    bool do_blocking;

    bool use_jastrow;
    bool use_coulomb;

    std::string system; //!< String specifing the system type, e.g. "Atoms".
    std::string sampling; //!< String specifing the sampling type, e.g. "IS".
    
    std::string runpath; //!< The directory which the simulation is set to run.

};

//! Struct used to initialize Minimization parameters.
/*!
 * \see ASGD 
 */
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
    
    arma::rowvec alpha; //!< Initial condition for the spatial variational parameter(s).
    arma::rowvec beta; //!< Initial condition for the Jastrow variational parameter(s).

};

//! Struct used to initialize output parameters.
struct OutputParams {
    bool dist_out; //!< If true, distributions are calculated for VMC/DMC.
    bool dmc_out; //!< If true, DMC outputs data to file each cycle.
    bool ASGD_out; //!< If true, ASGD outputs data to file each cycle.
};

//Forward declarations
class Orbitals;
class Potential;
class System;
class Sampling;
class Jastrow;

//! Struct used to initialize system objects.
/*!
 * The memory addresses allocated here will not change throughout the run.
 * \see Orbitals, Potential, System, Sampling, Jastrow
 */
struct SystemObjects {
    
    Orbitals* SP_basis;
    Potential* onebody_pot;
    System* SYSTEM;
    Sampling* sample_method;
    Jastrow* jastrow;

};


/*! \brief Class for handling standard output.
 * Only the master node has this object. 
 */
class STDOUT {
public:
    
    virtual void cout(std::stringstream & a) {
        std::cout << a.str() << std::endl;
        a.str(std::string());
        a.clear();
    }
};
 
/*! \brief Class for suppressing standard output.
 * Every node but the master has this. If-tests around cout is avoided.
 */
class NO_STDOUT : public STDOUT {
public:

    virtual void cout(std::stringstream & a) {
        a.str(std::string());
        a.clear();
    }
};

#include "Sampler/Sampler.h"
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
#include "Orbitals/DiTransform/DiTransform.h"
#include "Orbitals/ExpandedBasis/ExpandedBasis.h"

//#include "HartreeFock/HartreeFock.h"

#include "Potential/Potential.h"
#include "Potential/Coulomb/Coulomb.h"
#include "Potential/Harmonic_osc/Harmonic_osc.h"
#include "Potential/AtomCore/AtomCore.h"
#include "Potential/DiAtomCore/DiAtomCore.h"
#include "Potential/DoubleWell/DoubleWell.h"

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

