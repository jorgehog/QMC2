/* 
 * File:   ExpandedBasis.cpp
 * Author: jorgmeister
 * 
 * Created on October 9, 2012, 5:09 PM
 */

#include "../../QMCheaders.h"

ExpandedBasis::ExpandedBasis(GeneralParams & gp, Orbitals* basis, int basis_size, std::string coeffPath)
: Orbitals(gp.n_p, gp.dim) {

    this->basis = basis;
    this->basis_size = basis_size;
    coeffs = arma::zeros<arma::mat > (n_p, basis_size);

    std::ifstream coeffs_file;
    coeffs_file.open(coeffPath.c_str());

    for (int i = 0; i < n_p; i++) {
        for (int m = 0; m < basis_size; m++) {
            coeffs_file >> coeffs(i, m);
        }
    }

    coeffs_file.close();

}

double ExpandedBasis::phi(const Walker* walker, int particle, int q_num) const {

    double value = 0;
    for (int m = 0; m < basis_size; m++) {
        value += coeffs(particle, m) * basis->phi(walker, particle, q_num);
    }

    return value;

}

double ExpandedBasis::del_phi(const Walker* walker, int particle, int q_num, int d) const {

    double value = 0;
    for (int m = 0; m < basis_size; m++) {
        value += coeffs(particle, m) * basis->del_phi(walker, particle, q_num, d);
    }

    return value;

}

double ExpandedBasis::lapl_phi(const Walker* walker, int particle, int q_num) const {

    double value = 0;
    for (int m = 0; m < basis_size; m++) {
        value += coeffs(particle, m) * basis->lapl_phi(walker, particle, q_num);
    }

    return value;

}