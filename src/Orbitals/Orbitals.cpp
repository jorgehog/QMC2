/* 
 * File:   Orbitals.cpp
 * Author: jorgehog
 * 
 * Created on 26. juni 2012, 17:41
 */

#include "../QMCheaders.h"

Orbitals::Orbitals(int n_p, int dim) {
    this->n_p = n_p;
    this->n2 = n_p / 2;
    this->dim = dim;

    h = 1E-4;
    h2 = 1 / (h * h);
    two_h = 1 / (2 * h);
    

    max_implemented = 15; //for 30 particles
    basis_functions = new BasisFunctions*[max_implemented];

    dell_basis_functions = new BasisFunctions**[dim];
    for (int i = 0; i < dim; i++) {
        dell_basis_functions[i] = new BasisFunctions*[max_implemented];
    }

    lapl_basis_functions = new BasisFunctions*[max_implemented];
}

Orbitals::Orbitals() {

}

double Orbitals::num_diff(const Walker* walker, int particle, int q_num, int d) {

    Walker* diff_walker = new Walker(n_p, dim, false);
    diff_walker->r = walker->r;
    diff_walker->r(particle, d) += h;
    diff_walker->calc_r_i2(particle);
    set_qnum_indie_terms(diff_walker, particle);
    double phi_pluss = basis_functions[q_num]->eval(diff_walker, particle);

    diff_walker->r(particle, d) -= 2 * h;
    diff_walker->calc_r_i2(particle);
    set_qnum_indie_terms(diff_walker, particle);
    double phi_minus = basis_functions[q_num]->eval(diff_walker, particle);

    delete diff_walker;

    return (phi_pluss - phi_minus) / (2 * h);

}

double Orbitals::num_ddiff(const Walker* walker, int particle, int q_num) {
    Walker* diff_walker = new Walker(n_p, dim, false);
    diff_walker->r = walker->r;

    double ddiff = 0;
    for (int k = 0; k < dim; k++) {

        diff_walker->r(particle, k) += h;
        diff_walker->calc_r_i2(particle);
        set_qnum_indie_terms(diff_walker, particle);
        ddiff += basis_functions[q_num]->eval(diff_walker, particle);

        diff_walker->r(particle, k) -= 2 * h;
        diff_walker->calc_r_i2(particle);
        set_qnum_indie_terms(diff_walker, particle);
        ddiff += basis_functions[q_num]->eval(diff_walker, particle);

        diff_walker->r(particle, k) += h;
        diff_walker->r2(particle) = walker->r2(particle);
    }

    //contribution from phi(r, q_num) is constant in loop and already tabulated.
    ddiff = (ddiff - 2 * dim * walker->phi(particle, q_num)) * h2;

    delete diff_walker;

    return ddiff;
}

