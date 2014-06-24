#pragma once


#include <sstream>
#include <string>
#include <vector>
#include <armadillo>

#include "../Sampler/Sampler.h"


namespace QMC2
{


struct GeneralParams;
struct ParParams;
struct SystemObjects;

class ErrorEstimator;
class Walker;
class System;
class Jastrow;
class Sampling;
class Orbitals;
class Distribution;

/*! \brief The QMC superclass.
 * Holds implementations of general functions for both VMC and DMC in order to
 * avoid rewriting code and emphasize the similarities.
 */
class QMC {
protected:

    unsigned int output_tresh;

    unsigned int n_w_size; //!< The total number of allocated walkers.

    std::string runpath; //!< The directory which the simulation is set to run.

    std::string dist_path; //!< The path where the distribution are saved.
    arma::mat dist; //!< Matrix holding positional data for the distribution.
    int last_inserted; //!< Index of last inserted positional data.
    int dist_tresh; //!< Amount of cycles to skip in between storing position data.
    int dist_size;

    bool is_master;
    bool parallel;
    unsigned int node;
    unsigned int n_nodes;

    bool silent;

    unsigned int p_start; //!< Is always zero unless a deadlock is initialized.

    unsigned int n_c; //!< The number of Monte-Carlo cycles.
    unsigned int thermalization; //!< The number of thermalization steps.
    unsigned int cycle; //!< The current Monte-Carlo cycle.
    //! The number of walkers
    /*! VMC stores this many cycles in case of DMC */
    unsigned int n_w;

    unsigned int n_p;
    unsigned int n2;
    unsigned int dim;

    unsigned long int accepted; //!< Number of accepted moves.
    unsigned long int total_samples; //!< Total number of moves.

    double local_E; //!< The last calculated local energy.

    double dtOrig;

    double error;

    Walker *trial_walker; //!< The trial walker used to test a move.
    Walker **original_walkers; //!< A list of n_w walkers used in DMC.

    Jastrow *jastrow; //!< The Jastrow object.
    Sampling *sampling; //!< The Sampling object.
    System *system; //!< The system object.
    ErrorEstimator* error_estimator; //!< The error estimator.

    Distribution* distribution;

    std::vector<Sampler*> samplers;

    virtual void reset_all();

    void update_samplers(Walker *walker, double weight);

    //! Method for setting the trial position of the QMC method's walkers.
    /*!
     * \see Sampling::set_trial_pos()
     */
    virtual void set_trial_positions() = 0;

    //! Method for diffusing a walker one time step.
    /*!
     * The trial walker must equal the original walker in input.
     * The original walker is updated on output.
     */
    void diffuse_walker(Walker* original, Walker* trial);

    //! Method for calculating the acceptance ratio used in the Metropolis test.
    double get_acceptance_ratio(const Walker* walker_pre, const Walker* walker_post, int particle) const;

    //! Method for deciding whether or not to accept a move.
    /*!
     * Wraps the metropolis sampling with possibilities of overriding.
     * \see System::allow_transition()
     */
    virtual bool move_autherized(double A) = 0;

    //! Method for performing the metropolis test after when diffusing a walker.
    /*!
     * @param A The acceptance ratio calulated by get_acceptance_ratio().
     */
    bool metropolis_test(double A);

    //! Method for updating the walker after an accepted step.
    /*!
     * Given a particle number, the method only updates the changed values.
     * @param walker_post Walker at current time step
     * @param walker_pre Walker at previous time step
     */
    void update_walker(Walker*walker_pre, const Walker* walker_post, int particle) const;

    //! Method for reseting the walker after a rejected step.
    /*!
     * Given a particle number, the method only resets the changed values.
     * @param walker_post Walker at current time step
     * @param walker_pre Walker at previous time step
     */
    void reset_walker(const Walker* walker_pre, Walker* walker_post, int particle) const;

    //! Method for (hard) copying a walker object.
    /*!
     * @param parent,child The parent is copied to the child.
     */
    void copy_walker(const Walker* parent, Walker* child) const;

    //! Method for calculating the necessary quantities needed in order to calculate the local energy.
    /*!
     * \see Sampling::calculate_energy_necessities()
     */
    void calculate_energy_necessities(Walker* walker) const;

    //! Method for calculating the kinetic energy of a walker.
    double get_KE(const Walker* walker);

    void update_subsamples(double weight = 1.0);

    void push_subsamples();

    void finalize();

    //! Method for storing positional data.
    /*!
     * Stored in the dist matrix. Used by OutputHandler::Distribution.
     * \see Distribution::dump(), VMC::save_distribution(), DMC::save_distribution()
     */
    virtual void save_distribution() = 0;

    //! Method for performing node communication.
    virtual void node_comm() = 0;

    //!Method for calculating the distribution
    void finalize_distribution();
    
    /*!
     * Estimates and finalizes the ErrorEstimator object initialized in the error_estimator vector.
     */
    void estimate_error();

    /*!
     * Since the spatial wave function is split, certain values are unchanged if
     * the moved particle has opposite spin. Assuming a two-level system, the
     * first half of the particles are assumed to have one spin value, and the
     * second half the other.
     *
     * This method sets the start and end position of the block that needs to be
     * changed.
     *
     * \see System::start, System::end
     */
    void set_spin_state(unsigned int particle) const;

    
    //! Method used for testing the optimized ratio calculation.
    /*!
     * Compares to brute force computation of the wave function values.
     * @params R_qmc The optimized trial wave function ratio (spatial and Jastrow).
     */
    void test_ratios(const Walker* walker_pre, const Walker* walker_post, int particle, double R_qmc) const;

    //! Method for testing the optimized gradients calculation.
    /*!
     * Compares with finite difference calculation.
     */
    void test_gradients(Walker* walker);

    static std::string m_currentlyRunningMethod;
    
    
public:

    Sampler* kinetic_sampler;

    static void setCurrentlyRunningMethod(const std::string name)
    {
        m_currentlyRunningMethod = name;
    }

    static std::string currentlyRunningMethod()
    {
        return m_currentlyRunningMethod;
    }

    //! Constructor.
    /*!
     * @params K K times n_w walkers are initialized. K != 0 only sensible in DMC.
     */
    QMC(GeneralParams &, int n_c,
        SystemObjects &,
        ParParams &,
        bool silent,
        double dt,
        int n_w,
        int K = 1);

    //! Method used for executing the solver.
    virtual void run_method(bool initialize = true) = 0;

    void add_subsample(Sampler * sampler)
    {
        samplers.push_back(sampler);
    }

    void dump_subsamples(bool mean_of_means = false);

    virtual double get_energy() const = 0;

    double get_error() const
    {
        return error;
    }

    //! Method for calculating the Quantum Force.
    void get_QF(Walker* walker) const;

    //! Method for calculating the gradients after moving a particle.
    /*!
     * \see Jastrow::get_grad(), System::get_spatial_grad()
     */
    void get_gradients(const Walker* walker_pre, Walker* walker_post, int particle) const;

    //! Method for calculating the full gradients.
    /*!
     * \see Jastrow::get_grad(), System::get_spatial_grad()
     */
    void get_gradients(Walker* walker) const;

    //! Method for calculating the Laplacian of all walkers

    /*!
     * \see System::get_spatial_lapl_sum(), Jastrow::get_lapl_sum()
     */
    void get_laplsum(Walker* walker) const;

    //! Method for calculating the wave functions value at a given walker's position.

    /*!
     * \see System::get_spatial_wf(), Jastrow::get_val()
     */
    double get_wf_value(const Walker* walker) const;

    //! Method for calculating the local energy.

    /*!
     * \see get_KE(), System::get_potential_energy()
     */
    double calculate_local_energy(const Walker* walker);

    //! Method for calculating the acceptance ratio.
    void get_accepted_ratio();

    //! Cleans up initializations which distorts successive use of objects
    void clean();
    
    //! Method for setting the error estimator.
    void set_error_estimator(ErrorEstimator* error_estimator) {
        this->error_estimator = error_estimator;
    }

    //! Method for standard output.
    virtual void output() = 0;

    System* get_system_ptr() const {
        return system;
    }

    Sampling* get_sampling_ptr() const {
        return sampling;
    }

    Jastrow* get_jastrow_ptr() const {
        return jastrow;
    }

    Orbitals* get_orbitals_ptr() const;



};

}
