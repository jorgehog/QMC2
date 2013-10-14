#ifndef MOLECULARCOULOMB_H
#define MOLECULARCOULOMB_H

#include "../Potential.h"

struct GeneralParams;
class NBodyTransform;

class MolecularCoulomb : public Potential
{
public:
    MolecularCoulomb(GeneralParams & gP, NBodyTransform * orbitals);

    double get_pot_E(const Walker *walker) const;

private:

    NBodyTransform * nBodyOrbitals;

};

#endif // MOLECULARCOULOMB_H
