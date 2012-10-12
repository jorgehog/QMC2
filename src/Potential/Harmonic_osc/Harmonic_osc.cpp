/* 
 * File:   Harmonic_osc.cpp
 * Author: jorgmeister
 * 
 * Created on October 12, 2012, 2:42 PM
 */

#include "../../QMCheaders.h"


Harmonic_osc::Harmonic_osc(GeneralParams & gP)
: Potential(gP.n_p, gP.dim) {

    this->w = gP.w;

}


double Harmonic_osc::get_pot_E(const Walker* walker) const {

    double e_potential = 0;

    for (int i = 0; i < n_p; i++) {
        e_potential += 0.5 * w * w * walker->get_r_i2(i);
    }

    return e_potential;
}
