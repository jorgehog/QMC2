#include "oxygen3_21g.h"

#include "../BasisSetCodeMilad/splitValence/o_321g.h"


Oxygen3_21G::Oxygen3_21G() :
    GaussianFitted(8, 3)
{

    name = "O_321G";

    O_321G* basis = new O_321G();

    for (int i = 0; i < basis->getNumContracted(); ++i) {
        addGaussianFitFromCGTOs(basis->getContracted(i));
    }

}
