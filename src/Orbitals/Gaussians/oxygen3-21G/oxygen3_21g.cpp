#include "oxygen3_21g.h"

#include "../BasisSetCodeMilad/splitValence/o_321g.h"


Oxygen3_21G::Oxygen3_21G() :
    GaussianFitted(8, new O_321G)
{
    name = "O_321G";
}
