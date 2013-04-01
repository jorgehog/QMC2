/* 
 * File:   ASGD.h
 * Author: jorgehog
 *
 * Created on 23. August 2012, 16:52
 */

#ifndef ASGD_H
#define	ASGD_H

/*! \brief Implementation for the Adaptive Stochastic Gradient Descent method (ASGD)
 *  Used to find optimal variational parameters using adaptive step lengths.
 */
class ASGD : public Minimizer {
protected:

    int n_c; //!< The correlation length between storing two walkers after thermalization. 
    int n_c_SGD; //!< The number of samples used to estimate expectation values.
    int SGDsamples; //!< The number of ASGD cycles.
    int n_walkers; //!< The number of walkers.
    int thermalization; //!< The number of thermalization cycles used on walkers.
    
    int sample; //!< The current ASGD cycle.

    double t_prev; //!< The previous t.
    double t; //!< The current t.
    double step; //!< The current calculates step.
    double max_step; //!< The maximum threshold on a step.

    double E; //!< The energy summation variable used to calculate the mean.

    double a; //!< ASGD step parameter. 
    double A; //!< ASGD step parameter. 
    double f_min; //!< ASGD step parameter. 
    double f_max; //!< ASGD step parameter. 
    double w; //!< ASGD step parameter. 

    Walker** walkers; //!< The walkers used to sample expectation values.
    Walker** trial_walkers;

    arma::rowvec parameter; // UNUSED?

    arma::rowvec gradient; //!< Sumamtion vector for the trial wave function's variational derivatives.
    arma::rowvec gradient_local; //!< Sumamtion vector for the trial wave function's variational derivatives times the energy.

    arma::rowvec gradient_old; //!< The previous total gradient.
    arma::rowvec gradient_tot; //!< The current total gradient.
    
    //! Method for calculating the total gradient.
    /*!
     * Updates the error estimator with statistics.
     */
    void get_total_grad();
    
    //! Calculates the step and updates parameters.
    void update_parameters();
    
    //! Standard output of the progress.
    void output_cycle();
    
    //! Thermalizes a set of walkers before the main loop.
    void thermalize_walkers();

    //! Function for calculating the adaptive step.
    double f(double x) {
        return f_min + (f_max - f_min) / (1 - (f_max / f_min) * exp(-x / w));
    }
    
    //! Method for updating the vectors needed to calculate the total variational derivative.
    /*!
     * Calculates the single particle variational derivatives V and accumulates V and V*e_local.
     * @param e_local The local energy of the current walker at the current time step.
     */
    void get_variational_gradients(Walker* walker, double e_local);

public:
    ASGD(VMC*, MinimizerParams &, const ParParams &);

    void minimize();

    friend class stdoutASGD;
};


#endif	/* ASGD_H */