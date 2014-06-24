#include "Harmonic_osc.h"

#include "../../structs.h"
#include "../../Walker/Walker.h"


using namespace QMC2;

Harmonic_osc::Harmonic_osc(GeneralParams & gP)
: Potential(gP.n_p, gP.dim, "Oscillator") {

    this->w = gP.systemConstant;
    
}


double Harmonic_osc::get_pot_E(const Walker* walker) const {

    double e_potential = 0;

    for (int i = 0; i < n_p; i++) {
        e_potential += walker->get_r_i2(i);
    }

    return 0.5 * w * w * e_potential;
}
