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

    QMC* qmc;

    double h;
    double h2;
    double two_h;

    BasisFunctions** basis_functions;
    BasisFunctions*** dell_basis_functions;
    BasisFunctions** lapl_basis_functions;

    virtual double get_parameter(int n) = 0;
    virtual void set_parameter(double parameter, int n) = 0;

    virtual double get_variational_derivative(const Walker* walker, int n);

    double num_diff(const Walker* walker, int particle, int q_num, int d);
    double num_ddiff(const Walker* walker, int particle, int q_num);

    /*
     Test functions
     */
    
    void testLaplace(const Walker* walker, int particle, int q_num);
    void testDell(const Walker* walker, int particle, int q_num, int d);
    
    friend class Minimizer;
    friend class ASGD;
    friend class stdoutASGD;

public:

    Orbitals(int n_p, int dim);
    Orbitals();

    virtual void set_qnum_indie_terms(const Walker* walker, int i) {
    };

    virtual double phi(const Walker* walker, int particle, int q_num);

    virtual double del_phi(const Walker* walker, int particle, int q_num, int d);

    virtual double lapl_phi(const Walker* walker, int particle, int q_num);

    void set_qmc_ptr(QMC* qmc) {
        this->qmc = qmc;
    }

};


#endif	/* ORBITALS_H */
