#pragma once


#include "../QMC/VMC/VMC.h"

#include <vector>
#include <armadillo>
#include <string>


namespace QMC2
{


class Jastrow;
class Orbitals;


/*!\brief Class for minimization methods used to obtain optimal variational parameters.
 */
class Minimizer {
protected:
    
    int n_nodes;
    bool is_master;

    VMC* vmc; //!< Uses VMC methods to calculate stochastic variational gradients.

    std::stringstream s; 
    
    int Nspatial_params; //!< The number of variational parameters in the spatial trial wave function.
    int Njastrow_params; //!< The number of variational parameters in the Jastrow factor.
    int Nparams; //!< The total number of variational parameters.

    //! Method for updating the variational parameters based on the previous step.
    /*!
     * Needs to be implemented by a subclass.
     */
    virtual void update_parameters() = 0;

    void initializeParameters();

    arma::rowvec alpha0;
    arma::rowvec beta0;

public:

    //! Constructor. 
    /*!
     * @param vmc The VMC object used for storing variational parameters and calculating stochastic gradients.
     * @param alpha Vector of initial conditions of spatial variational parameters
     * @param beta Vector of initial conditions of Jastrow variational parameters
     */
    Minimizer(VMC* vmc, const ParParams &, const arma::rowvec & alpha0, const arma::rowvec & beta0);


    Orbitals* get_orbitals() {
        return vmc->get_orbitals_ptr();
    }

    Jastrow* get_jastrow() {
        return vmc->get_jastrow_ptr();
    }

    //! Method for executing the minimization main solver.
    virtual void minimize(bool initialize = true) = 0;

};

}
