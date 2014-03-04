#pragma once


#include "../Potential.h"


namespace QMC2
{


struct GeneralParams;

/*! \brief Implementation of the Coulomb potential. 1/r_{ij}
 */
class Coulomb : public Potential {
public:

    Coulomb(GeneralParams &);

    double get_pot_E(const Walker* walker) const;

};

}
