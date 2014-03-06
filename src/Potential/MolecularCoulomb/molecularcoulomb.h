#pragma once


#include "../Potential.h"


namespace QMC2
{


struct GeneralParams;
class NBodyTransform;

class MolecularCoulomb : public Potential
{
public:
    MolecularCoulomb(GeneralParams & gP);

    double get_pot_E(const Walker *walker) const;

private:

    NBodyTransform * nBodyOrbitals;

};

}
