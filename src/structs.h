#ifndef MISC_H
#define MISC_H

#include "defines.h"

#include <armadillo>
#include <iostream>

#include <sys/time.h>

//! Struct used to initialize DMC parameters.
struct DMCparams {
    int n_c = 1000; //!< The number of cycles
    int therm = 1000; //!< Thermalization cycles
    int n_b = 100; //!< Number of block samples pr. walker pr. cycle
    int n_w = 1000; //!< Number of walkers.
    double dt = 0.001; //!< Time step
};

//! Struct used to initialize VMC parameters.
struct VMCparams {
    int n_c = 1000000; //!< The number of cycles.
    double dt = 0.005; //!< The time step.
};

//! Struct used to initialize the varational parameters.
struct VariationalParams {
    double alpha = 0.987; //!< The spatial variational parameter.
    double beta = 0.398; //!< The Jastrow variational parameter.
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
    int n_p = 2; //!< The number of particles.
    int dim = 2; //!< The dimension.
    seed_type random_seed = -time(NULL); //!< The random number generator's seed.
    
    //! The constant used in systems
    /*!
     * e.g. charge for atoms and oscillator frequency for quantum dots.
     */
    double systemConstant = 1;

    double R = 1.4; //! Center of mass coordinate for diatmic systems.

    bool deadlock = false; //! If true, freezes one particle;
    double deadlock_x = 0; //! Position of the locked particle. y=z=0;
    
    std::string runpath = "../scratch/QMC_SCRATCH/"; //!< The directory which the simulation is set to run.

};

//! Struct used to initialize Minimization parameters.
/*!
 * \see ASGD 
 */
struct MinimizerParams {
    double max_step = 0.1;
    double f_max = 1.0;
    double f_min = -0.5;
    double omega = 0.8;
    double A = 60;
    double a = 0.5;
    
    int SGDsamples = 2000;
    int n_w = 10;
    int therm = 10000;
    int n_c = 1000;
    int n_c_SGD = 400;
    
    arma::rowvec alpha = arma::zeros(1, 1) + 0.5; //!< Initial condition for the spatial variational parameter(s).
    arma::rowvec beta = arma::zeros(1, 1) + 0.5; //!< Initial condition for the Jastrow variational parameter(s).

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
    System* system;
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

inline void initMPI(struct ParParams & parParams, int argc, char ** argv){
#ifdef MPI_ON
    int node, n_nodes;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &n_nodes);
    MPI_Comm_rank(MPI_COMM_WORLD, &node);

    parParams.n_nodes = n_nodes;
    parParams.node = node;
    parParams.parallel = (parParams.n_nodes > 1);
    parParams.is_master = (parParams.node == 0);
#else
    parParams.parallel = false;
    parParams.node = 0;
    parParams.n_nodes = 1;
    parParams.is_master = true;
#endif
}

inline void scaleWithProcs(struct ParParams & parParams,
                           struct GeneralParams & generalParams,
                           struct MinimizerParams & minimizerParams,
                           struct VMCparams & vmcParams,
                           struct DMCparams & dmcParams){

    generalParams.random_seed -= parParams.node;
    minimizerParams.n_c_SGD /= parParams.n_nodes;

    vmcParams.n_c /= parParams.n_nodes;
    dmcParams.n_w /= parParams.n_nodes;

}







#endif /* MISC_H */
