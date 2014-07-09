#pragma once

#include <functional>
#include <vector>
#include <string>


namespace QMC2
{


class Walker;
class Orbitals;
class Potential;

/*! \brief The system class separating Fermions and Bosons.
 * Designed to generalize the solver in terms of particle species.
 */
class System {
protected:
    int n_p;
    int n2;
    int dim;

    //! The start point of separable calculations.
    /*!
     * Either 0 or N/2.
     * Since the spatial wave function is split, particles with spin
     * not equal that of the moved particle is unchanged and does not
     * need to be recalculated.
     */
    int start;

    //! The end point of separable calculations.
    /*!
     * Either N/2 or N.
     * Since the spatial wave function is split, particles with spin
     * not equal that of the moved particle is unchanged and does not
     * need to be recalculated.
     */
    int end;

    std::vector<Potential*> potentials; //!< A vector of potentials.
    Orbitals *orbital; //!< The single particle wave function object.

public:
    System(int n_p, int dim, Orbitals* orbital);
    System();

    //! Initializes the system before the main solver loop starts.
    /*!
     * Called by the Sampling class when trial positions are set.
     */
    virtual void initialize(Walker* walker) = 0;

    //! Method for adding a potential to the system.
    void add_potential(Potential* pot);

    //! Method for calculating the total potential energy.
    /*!
     * Iterates over all objects in the potentials vector and accumulates
     * their potential energies for the given walker.
     */
    double get_potential_energy(const Walker* walker);

    // Method for updating the system specific parts of a walker.
    /*!
     * \see QMC::update_walker()
     */
    virtual void update_walker(Walker* walker_pre, const Walker* walker_post, int particle) const = 0;

    // Method for resetting the system specific parts of a walker.
    /*!
     * \see QMC::reset_walker()
     */
    virtual void reset_walker(const Walker* walker_pre, Walker* walker_post, int particle) const = 0;

    //! Method for calculating the necessary values needed by the walker after a new step is made.
    /*!
     * Given a particle number, the method does not recompute unchanged values.
     * @param walker_old Walker at current time step.
     * @param walker_new Walker at previous time step.
     */
    virtual void calc_for_newpos(const Walker* walker_old, Walker* walker_new, int particle) = 0;

    //! Method for calculating the spatial wave function ratios between two subsequent time steps.
    virtual double get_spatial_ratio(const Walker* walker_pre, const Walker* walker_post, int particle) = 0;

    //! Method for calculating the spatial wave function's value at a given walkers position.
    virtual double get_spatial_wf(const Walker* walker) = 0;

    //! Method for calculating the changed part of the spatial gradient.
    /*!
     * Depending on which particle we moved, one of the spatial wave function parts (it is split)
     * will be unchanged.
     */
    virtual void get_spatial_grad(Walker* walker, int particle) const = 0;

    //! Method for calculating the full spatial gradient.
    virtual void get_spatial_grad_full(Walker* walker) const = 0;

    //! Method for calculating the sum of all Laplacians for a given walker.
    virtual double get_spatial_lapl_sum(Walker* walker) const = 0;

    //! Method for copying the system specific parts of a walker.
    /*!
     * \see QMC::copy_walker()
     */
    virtual void copy_walker(const Walker* parent, Walker* child) const = 0;

    //! Method allowing the system to override the Metropolis test.
    /*
     * Depending on the system at hand, certain scenarios may result in
     * absolute rejection of the step (fixed node approximation for Fermions)
     */
    virtual bool allow_transition() = 0;

    void forEachPotentialDo(std::function<void(Potential *pot)> func);

    void update_potential_samples(double weight = 1.0);
    
    void push_potential_samples();

    std::string dump_samples();

    void initializeSamplingErrorEstimators(int type, int n_c);

    void reset_potential_samples();
    
    void finalize_potential_samples();

    Orbitals* get_orbital_ptr() {
        return orbital;
    }

    /*!
     * \see QMC::set_spin_state()
     */
    void set_spin_state(int start, int end) {
        this->start = start;
        this->end = end;
    }
};

}
