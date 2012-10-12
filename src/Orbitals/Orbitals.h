/* 
 * File:   Orbitals.h
 * Author: jorgehog
 *
 * Created on 26. juni 2012, 17:41
 */

#ifndef ORBITALS_H
#define	ORBITALS_H

class Orbitals {
protected:
    int n_p;
    int n2;
    int dim;
    int max_implemented;

    BasisFunctions** basis_functions;
    BasisFunctions*** dell_basis_functions;
    BasisFunctions** lapl_basis_functions;

    virtual double get_parameter(int n) = 0;
    virtual void set_parameter(double parameter, int n) = 0;
    virtual double get_variational_derivative(const Walker* walker, int n) const = 0;

    friend class Minimizer;
    friend class ASGD;

public:

    Orbitals(int n_p, int dim);
    Orbitals();

    virtual double phi(const Walker* walker, int particle, int q_num) const;
    virtual double del_phi(const Walker* walker, int particle, int q_num, int d) const;
    virtual double lapl_phi(const Walker* walker, int particle, int q_num) const;

};


#endif	/* ORBITALS_H */
