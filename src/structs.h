#ifndef MISC_H
#define MISC_H

#include "defines.h"

#include <armadillo>
#include <iostream>

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

    bool deadlock; //! If true, freezes one particle;
    double deadlock_x; //! Position of the locked particle. y=z=0;
    
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





#endif /* MISC_H */