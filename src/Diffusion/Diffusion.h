/* 
 * File:   Diffusion.h
 * Author: jorgehog
 *
 * Created on 16. april 2012, 14:03
 */

#ifndef DIFFUSION_H
#define	DIFFUSION_H

#include "../defines.h"
class QMC;
class Walker;

/*! \brief Class containing rules for walker movement based on diffusion models.
 * Serves as class member in the Sampling class. Brute force implies the Simple
 * diffusion model, while Importance Sampling implies the Fokker Planck diffusion.
 * @see Brute_Force, Importance.
 */
class Diffusion {
protected:
    int n_p;
    int dim;

    double timestep; //!< The discrete time step
    double D; //!< The diffusion constant
    long random_seed; //!< The random seed. Needs to be stored for some RNGs to work.

    double std; //!< The standard deviation from QMC stored for efficiency. sqrt(2D*timestep).


public:

    Diffusion(int n_p, int dim, double timestep, seed_type random_seed, double D);

    //! Virtual function returning the new position.
    /*!
     * Returns the simple diffusion step if not overridden.
     * @param i Particle number.
     * @param j dimension (x,y,z).
     * @return The new position (relative to the old).
     */
    virtual double get_new_pos(const Walker* walker, int i, int j);

    //! Calculates the Diffusion Green's function ratio needed by metropolis.
    /*!
     * @param walker_post Walker at current time step.
     * @param walker_pre Walker at previous time step.
     * @return The Diffusion Green's function ratio.
     */
    virtual double get_g_ratio(const Walker* walker_post, const Walker* walker_pre) const = 0;

    //! Calls a uniform random number generator.
    /*!
     * Returns a random uniform number on [0,1).
     */
    double call_RNG();

    //! Function for altering the time step. 
    /*!
     * Takes care of consequences. Time step should only be altered using this
     * function.
     */
    void set_dt(double dt);

    double get_dt() const {
        return timestep;
    }

    double get_std() const {
        return std;
    }

};


#endif	/* DIFFUSION_H */

