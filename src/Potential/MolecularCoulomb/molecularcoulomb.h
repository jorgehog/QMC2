#pragma once


#include "../Potential.h"

#include <vector>
#include <armadillo>

using namespace std;
using namespace arma;

namespace QMC2
{


struct GeneralParams;

class MolecularCoulomb : public Potential
{
public:
    MolecularCoulomb(GeneralParams & gP, const vector<const rowvec*> &corePositions, const rowvec &coreCharge);

    double get_pot_E(const Walker *walker) const;

    void makeRRelNucleiMatrix();

    const uint NCores;
    const vector<const rowvec*> corePositions;
    const rowvec coreCharge;

    mat r_rel_nuclei;

};

}
