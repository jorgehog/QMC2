/* 
 * File:   Potential.cpp
 * Author: jorgehog
 * 
 * Created on 13. april 2012, 17:21
 */

#include "../QMCheaders.h"

Potential::Potential(int n_p, int dim) {

    this->n_p = n_p;
    this->dim = dim;

}

Harmonic_osc::Harmonic_osc(int n_p, int dim, double w)
: Potential(n_p, dim) {

    this->w = w;

}

double Harmonic_osc::get_pot_E(const Walker* walker) const {

    double e_potential = 0;

    for (int i = 0; i < n_p; i++) {
        e_potential += 0.5 * w * w * walker->get_r_i2(i);
    }

    return e_potential;
}

Coulomb::Coulomb(int n_p, int dim)
: Potential(n_p, dim) {

}

double Coulomb::get_pot_E(const Walker* walker) const {
    
    double e_coulomb = 0;

    for (int i = 0; i < n_p - 1; i++) {
        for (int j = i + 1; j < n_p; j++) {
            e_coulomb += 1 / walker->r_rel(i, j);
        }
    }

    return e_coulomb;
}
