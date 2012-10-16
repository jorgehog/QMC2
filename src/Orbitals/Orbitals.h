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

    double h;
    double h2;
    double two_h;

    BasisFunctions** basis_functions;
    BasisFunctions*** dell_basis_functions;
    BasisFunctions** lapl_basis_functions;

    virtual double get_parameter(int n) = 0;
    virtual void set_parameter(double parameter, int n) = 0;
    virtual double get_variational_derivative(const Walker* walker, int n) const = 0;

    double num_diff(const Walker* walker, int particle, int q_num, int d) const;
    double num_ddiff(const Walker* walker, int particle, int q_num) const;

    friend class Minimizer;
    friend class ASGD;

public:

    Orbitals(int n_p, int dim);
    Orbitals();

    virtual double phi(const Walker* walker, int particle, int q_num) const {
        return basis_functions[q_num]->eval(walker, particle);
    }

    virtual double del_phi(const Walker* walker, int particle, int q_num, int d) const {
        return dell_basis_functions[d][q_num]->eval(walker, particle);
    }

    virtual double lapl_phi(const Walker* walker, int particle, int q_num) const {
        return lapl_basis_functions[q_num]->eval(walker, particle);
    }

};


#endif	/* ORBITALS_H */
