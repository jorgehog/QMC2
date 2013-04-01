/* 
 * File:   MINIMIZER.h
 * Author: jorgehog
 *
 * Created on 23. August 2012, 16:52
 */

#ifndef MINIMIZER_H
#define	MINIMIZER_H

/*!\brief Class for minimization methods used to obtain optimal variational parameters.
 */

class Minimizer {
protected:
    
    int n_nodes;
    bool is_master;

    VMC* vmc; //!< Uses VMC methods to calculate stochastic variational gradients.

    STDOUT* std_out; //!< Output object. Wraps and replaces std::cout.
    std::stringstream s; 
    
    int Nspatial_params; //!< The number of variational parameters in the spatial trial wave function.
    int Njastrow_params; //!< The number of variational parameters in the Jastrow factor.
    int Nparams; //!< The total number of variational parameters.

    std::vector<OutputHandler*> output_handler; //!< Either contains a stdoutASGD object or not.
    std::vector<ErrorEstimator*> error_estimators; //!< One ErrorEstimator object pr. variational parameter. 

    /*!
     * Iterates over the output objects in the output_handler vector. No if-tests.
     */
    void dump_output();
    
    /*!
     * Calls the finalize function for the object in the output_handler vector.
     */
    void finalize_output();
    
    /*!
     * Estimates and finalizes the ErrorEstimator objects initialized in the error_estimators vector.
     */
    void error_output();

    //! Method for updating the variational parameters based on the previous step.
    /*!
     * Needs to be implemented by a subclass.
     */
    virtual void update_parameters() = 0;

    //virtual void get_variational_gradients(); is this a general req?

public:

    //! Constructor. 
    /*!
     * @param vmc The VMC object used for storing variational parameters and calculating stochastic gradients.
     * @param alpha Vector of initial conditions of spatial variational parameters
     * @param beta Vector of initial conditions of Jastrow variational parameters
     */
    Minimizer(VMC* vmc, const ParParams &, const arma::rowvec & alpha, const arma::rowvec & beta);

    //! Method used for loading the stdoutASGD object.
    void add_output(OutputHandler* output_handler);

    Orbitals* get_orbitals() {
        return vmc->get_orbitals_ptr();
    }

    Jastrow* get_jastrow() {
        return vmc->get_jastrow_ptr();
    }

    //! Method for executing the minimization main solver.
    virtual void minimize() = 0;

    //! Method for dumping variational parameter values to screen.
    void output(std::string message, double number = -1);

    //! Method used to add error estimators.
    void add_error_estimator(ErrorEstimator* error_estimator) {
        this->error_estimators.push_back(error_estimator);
    }
};


#endif	/* MINIMIZER_H */