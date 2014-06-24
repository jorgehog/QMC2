/* 
 * File:   Coulomb.cpp
 * Author: jorgmeister
 * 
 * Created on October 12, 2012, 2:41 PM
 */

#include "Coulomb.h"

#include "../../structs.h"
#include "../../Walker/Walker.h"


using namespace QMC2;


Coulomb::Coulomb(GeneralParams & gP)
: Potential(gP.n_p, gP.dim, "Coulomb") {

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

