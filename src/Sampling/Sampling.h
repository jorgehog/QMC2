#pragma once


#include "../Diffusion/Diffusion.h"
#include "../QMC/QMC.h"


namespace QMC2
{


class Sampling {
protected:
    int n_p;
    int n2;
    int dim;

    /*!
     * \see System::start
     */
    int start;

    /*!
     * \see System::end
     */
    int end;

    bool deadlock; //! If true, a particle is frozen in the simulation.
    double deadlock_x; //! Position of the frozen particle. y=z=0.

    void set_deadlock(const double deadlock_x);
    void clear_deadlock();

    //!The Diffusion object.
    /*! \see Diffusion */
    Diffusion* diffusion;
    QMC* qmc; //!< The QMC main solver object. Needed to access e.g. the system object.

    friend QMC::QMC(GeneralParams &, int,
                    SystemObjects &,
                    ParParams &,
                    bool,
                    double,
                    int,
                    int);
    
    friend void QMC::clean();

public:

    Sampling(int n_p, int dim);
    Sampling();

    //! Method for updating the position of a walker's particle.
    /*!
     * Sets a new position according to the diffusion rules,
     * and calls all the functions necessary to get all the values updates,
     * e.g. System::calc_for_new_pos()
     */
    void update_pos(const Walker* walker_pre, Walker* walker_post, int particle) const;

    //! Method for updating the sampling specific necessary values.
    virtual void update_necessities(const Walker* walker_pre, Walker* walker_post, int particle) const = 0;

    // Method for updating the sampling specific parts of a walker.
    /*!
     * \see QMC::update_walker()
     */
    virtual void update_walker(Walker* walker_pre, const Walker* walker_post, int particle) const = 0;

    // Method for resetting the sampling specific parts of a walker.
    /*!
     * \see QMC::reset_walker()
     */
    virtual void reset_walker(const Walker* walker_pre, Walker* walker_post, int particle) const = 0;

    //! Method for setting the trial position for a given walker.
    void set_trial_pos(Walker* walker);

    //! Method for setting up the single particle orbitals and it's derivatives for a given walker.
    void set_trial_states(Walker* walker); //## Should be Fermion specific.

    //! Method for calculating the sampling specific necessary values.
    /*!
     * Called after a trial position is set.
     */
    virtual void get_necessities(Walker* walker) = 0;

    //! Method for calculating the sampling specific necessary values in order to compute the local energy.
    virtual void calculate_energy_necessities(Walker* walker) const = 0;

    //! Method for copying the sampling specific parts of a walker.
    /*!
     * \see QMC::copy_walker()
     */
    virtual void copy_walker(const Walker* parent, Walker* child) const = 0;

    //! Method for calculating the diffusion Green's function ratios.

    /*!
     * See the Diffusion class for documentation.
     */
    virtual double get_g_ratio(const Walker* walker_post, const Walker* walker_pre) const {
        return diffusion->get_g_ratio(walker_post, walker_pre);
    }

    //! Calculates the Branching Green's function ratio needed by DMC.

    /*! 
     * @param E_x Energy at current time step
     * @param E_y Energy at previous time step
     * @return The Branching Green's function ratio
     */
    double get_branching_Gfunc(double E_x, double E_y, double E_T) const {
        return exp(-(0.5 * (E_x + E_y) - E_T) * diffusion->get_dt());
    }

    void set_qmc_ptr(QMC* qmc) {
        this->qmc = qmc;
    }

    void set_dt(double dt) {
        diffusion->set_dt(dt);
    }

    double get_dt() const {
        return diffusion->get_dt();
    }

    double get_std() const {
        return diffusion->get_std();
    }

    //! Calls a uniform random number generator.

    /*!
     * Returns a random uniform number on [0,1).
     */
    double call_RNG() {
        return diffusion->call_RNG();
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
