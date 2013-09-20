/* 
 * File:   ExpandedBasis.cpp
 * Author: jorgmeister
 * 
 * Created on October 9, 2012, 5:09 PM
 */

#include "ExpandedBasis.h"

#include "../../structs.h"
#include "../../Walker/Walker.h"


ExpandedBasis::ExpandedBasis(GeneralParams & gp, Orbitals* basis, int basis_size)
: Orbitals(gp.n_p, gp.dim) {

    name = "Expanded" + basis->getName();

    this->basis = basis;
    this->basis_size = basis_size;
    coeffs = arma::zeros<arma::mat > (n_p, basis_size);

    //load coeffs in some intelligent way.

}

double ExpandedBasis::phi(const Walker* walker, int particle, int q_num) {

    double value = 0;

    //Dividing basis_size by half assuming a two-level system.
    //In case of Bosons, expanding s.p. w.f. does not make sence.
    for (int m = 0; m < basis_size / 2; m++) {
        value += coeffs(q_num, m) * basis->phi(walker, particle, m);
    }

    return value;

}

double ExpandedBasis::del_phi(const Walker* walker, int particle, int q_num, int d) {

    double value = 0;
    for (int m = 0; m < basis_size; m++) {
        value += coeffs(particle, m) * basis->del_phi(walker, particle, q_num, d);
    }

    return value;

}

double ExpandedBasis::lapl_phi(const Walker* walker, int particle, int q_num) {

    double value = 0;
    for (int m = 0; m < basis_size; m++) {
        value += coeffs(particle, m) * basis->lapl_phi(walker, particle, q_num);
    }

    return value;

}
