/* 
 * File:   Orbitals.cpp
 * Author: jorgehog
 * 
 * Created on 26. juni 2012, 17:41
 */

#include "../QMCheaders.h"

Orbitals::Orbitals(int n_p, int dim) {
    this->n_p = n_p;
    this->n2 = n_p/2;
    this->dim = dim;

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

double Orbitals::phi(const Walker* walker, int particle, int q_num) const {
    return basis_functions[q_num]->eval(walker, particle);
}

double Orbitals::del_phi(const Walker* walker, int particle, int q_num, int d) const {
    return dell_basis_functions[d][q_num]->eval(walker, particle);
}

double Orbitals::lapl_phi(const Walker* walker, int particle, int q_num) const {
    return lapl_basis_functions[q_num]->eval(walker, particle);
}


