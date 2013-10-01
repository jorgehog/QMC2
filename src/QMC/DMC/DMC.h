/* 
 * File:   DMC.h
 * Author: jorgmeister
 *
 * Created on October 12, 2012, 2:42 PM
 */

#ifndef DMC_H
#define	DMC_H

#include "../QMC.h"
class VMC;
struct DMCparams;
class stdoutDMC;

/*! \brief Implementation of the Diffusion Monte-Carlo Method.
 * Very little needs to be added when the QMC superclass holds all the
 * general functionality.
 */
class DMC : public QMC {
protected:

    bool thermalized; //!< Flag used to indicate whether to start sampling or not.
    bool initializedFromVMC;

    unsigned int n_w_last; //!< The amount of walkers at the time the walker loop was initiated.
    unsigned int n_w_tot; //!< The total number of walkers across all nodes.
    arma::uvec n_w_list; //!< List of the number of walkers of each node.

    bool force_comm; //!< Flag set true if population should be renormalized.
    
    unsigned int deaths; //< The number of walkers who died the last cycle.

    unsigned int block_size; //< The number of block samples for each walker for each cycle.
    unsigned int samples; //< The number of samples to the expectation value made every cycle.

    double dmc_E; //!< The DMC energy.
    double dmc_E_unscaled; //!< The accumulative DMC energy: The sum of all previous trial energies.

    double E_T; //!< The trial energy at the current cycle.
    double E; //!< The accumulative energy for each cycle.

    stdoutDMC* DMCout;


    //! Method for setting the trial position of all DMC walkers.
    /*!
     * In case VMC is not run prior to DMC, trial positions must be set. 
     */
    void set_trial_positions();

    //! Method for iterating a walker one time step.
    /*!
     * @param thermalized Flag to indicate whether to start sampling or not.
     * @param k The index of the walker.
     */
    void iterate_walker(int k);

    //! Method for the birth/death process of a walker.
    /*!
     * @param GB The branching Green's Function.
     * @param k The index of the walker.
     */
    void Evolve_walker(int k, double GB);

    //! Method for cleaning up the dead walkers and compact the list.
    void bury_the_dead();

    //! Method for updating the DMC energy and calculating the new trial energy.
    void update_energies();

    /*!
     * In case of DMC, we must let the system have the possibility to
     * override the metropolis test (fixed node approximation in case of a
     * Fermion system) 
     */
    bool move_autherized(double A);

    void reset_parameters() {
        n_w_last = n_w;
        E = samples = deaths = 0;
    }

    /*!
     * For each process:
     * -Calculates the total number of walkers.
     * -Sums up the energies sampled.
     * -Sums up the total number of samples made. 
     */
    void node_comm();

    //! Method for storing positional data for the Distribtuon.
    /*!
     * Stores the position data from all currently alive DMC walkers 
     * every dist_tresh cycle.
     */
    void save_distribution();

    //! Method for sending a walker between two nodes.
    /*!
     * @param root The node from which the walker is sent.
     * @param root_id The index of the walker being sent from root.
     * @param dest The node which receives the walker.
     * @param dest_id The index where the walker is to be received.
     * \see Walker::send_soul(), Walker::recv_soul()
     */
    void switch_souls(unsigned int root, unsigned int root_id, unsigned int dest, unsigned int dest_id);

    //! Method for evening out the number of walkers on each node.
    void normalize_population();

    //! Method which deletes all walkers.
    void free_walkers();

    void reset_all();

public:

    //! Factor of empty space for walkers over initial walkers.
    /*! \see QMC::QMC() */
    static const unsigned int K = 50;
    static const unsigned int check_thresh = 25; //< The amount of cycles between every time the population is normalized.
    static const unsigned int sendcount_thresh = 20; //< Minimum threshold for initializing a population normalization.

    //! Constructor.
    DMC(GeneralParams &, DMCparams &, SystemObjects &, ParParams &, bool silent = false);

    void run_method(bool initialize = true);

    void initFromVMC(VMC* vmc);

    void output();

    friend class stdoutDMC;

};

#endif	/* DMC_H */
