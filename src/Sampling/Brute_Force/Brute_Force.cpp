#include "Brute_Force.h"

#include "../../structs.h"

#include "../../Diffusion/Simple/Simple.h"


using namespace QMC2;

Brute_Force::Brute_Force(GeneralParams & gP)
: Sampling(gP.n_p, gP.dim) {
    diffusion = new Simple(n_p, dim, 0, gP.random_seed);

}
