/* 
 * File:   Orbitals.h
 * Author: jorgehog
 *
 * Created on 26. juni 2012, 17:41
 */

#ifndef ORBITALS_H
#define	ORBITALS_H

#include "../HartreeFock/HartreeFock.h"

/*! \brief Superclass for the single particle orbital classes.
 * Handles everything specific regarding choice of single particle basis.
 */
class Orbitals {
protected:
    int n_p;
    int n2;
    int dim;

    int max_implemented; //!< The maximum number basis size supported for any system ##RYDD OPP

    QMC* qmc; //!< A pointer to the QMC solver object. Needed for numerical variational derivatives. 

    double h; //!< The step length for finite difference schemes.
    double h2;
    double two_h;

    arma::imat qnums; //!< Quantum number matrix needed by Hartree-Fock and the variational derivatives.

    BasisFunctions** basis_functions; //!< A vector maping a quantum number index to a single particle wave function.
    BasisFunctions*** dell_basis_functions; //!< A maxtrix maping a quantum number- and dimension index to a single particle wave function derivative.
    BasisFunctions** lapl_basis_functions; //!< A vector maping a quantum number index to a single particle wave function Laplacian.

    //! A method for retrieving variational parameters.
    /*! 
     * @param n Index of the sought variational parameter.
     */
    virtual double get_parameter(int n) = 0;

    //! A method for setting variational parameters.
    /*!
     * @param n Index of the sought variational parameter.
     * @param parameter The new value of the variational parameter.
     */
    virtual void set_parameter(double parameter, int n) = 0;

    //! A method for calculating the variational derivative. 
    /*!
     * By default uses a finite difference scheme.
     * Can be overridden to evaluate a closed form expression.
     * @param n Index of the sought variational parameter.
     */
    virtual double get_variational_derivative(const Walker* walker, int n);

    //! Method for calculating the single particle derivative using a finite difference scheme.
    /*!
     * Method del_phi() can be overridden to use this method in case no closed form expressions are implemented.
     * @param d The dimension for which the derivative should be calculated (x,y,z).
     */
    double num_diff(const Walker* walker, int particle, int q_num, int d);

    //! Method for calculating the single particle Laplacian using a finite difference scheme.
    /*!
     * Method lapl_phi() can be overridden to use this method in case no closed form expressions are implemented.
     * @param q_num The quantum number index.
     */
    double num_ddiff(const Walker* walker, int particle, int q_num);

    //! Method for validating closed form expressions for laplacians by comparing them to numerical calculations.
    /*!
     * @param q_num The quantum number index.
     */
    void testLaplace(const Walker* walker, int particle, int q_num);

    //! Method for validating closed form expressions for derivatives by comparing them to numerical calculations.
    /*!
     * @param q_num The quantum number index.
     * @param d The dimension for which the derivative should be calculated (x,y,z).
     */
    void testDell(const Walker* walker, int particle, int q_num, int d);

    //! Method for calculating the anti-symmetrized Coulumb matrix elements
    /*!
     * Used by Hartree Fock
     */

    virtual double get_coulomb_element(const arma::uvec & qnum_set);
    virtual double get_sp_energy(int qnum) const;
    


public:

    Orbitals(int n_p, int dim);
    Orbitals();

    //! Calculates single particle wave function terms which are independent of the quantum numbers
    /*!
     * If a term in the single particle functions are independent of the quantum number,
     * this function can be overridden to calculate them beforehand (for each particle), 
     * and rather extract the value instead of recalculating.
     * @param i Particle number.
     */
    virtual void set_qnum_indie_terms(Walker* walker, int i) {
    };

    //! Calculates the single particle wave function for a given walker's particle.
    /*!
     * @param q_num The quantum number index.
     */
    virtual double phi(const Walker* walker, int particle, int q_num);

    //! Calculates the single particle wave function derivative for a given walker's particle and dimension.
    /*!
     * @param q_num The quantum number index.
     * @param d The dimension for which the derivative should be calculated (x,y,z).
     */
    virtual double del_phi(const Walker* walker, int particle, int q_num, int d);

    //! Calculates the single particle wave function for a given walker's particle.
    /*!
     * @param q_num The quantum number index.
     */
    virtual double lapl_phi(const Walker* walker, int particle, int q_num);

    void set_qmc_ptr(QMC* qmc) {
        this->qmc = qmc;
    }
    
    friend class HartreeFock;
    friend class Minimizer;
    friend class ASGD;
    friend class stdoutASGD;
};


#endif	/* ORBITALS_H */
