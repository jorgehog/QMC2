#include "Diffusion.h"

#include "../QMC/QMC.h"

#ifdef RNG_ZIG
#include "RNGs/zignor.h"
#include "RNGs/zigrandom.h"
#endif

#ifdef RNG_NUMREC
#include "RNGs/gaussian_deviate.h"
#include "RNGs/ran2.h"
#endif


using namespace QMC2;

Diffusion::Diffusion(int n_p, int dim, double timestep, seed_type random_seed, double D) {
    this->n_p = n_p;
    this->dim = dim;
    this->timestep = timestep;
    this->random_seed = random_seed;
    this->D = D;
    this->std = sqrt(2 * D * timestep);


#ifdef RNG_ZIG
    int cseed = 100;
    int random_seed2 = random_seed * 3;
    RanSetSeed_MWC8222(&random_seed2, cseed);
    RanNormalSetSeedZig32(&random_seed, 5);
#endif

#ifdef RNG_NUMREC
    if (random_seed < 0) {
        random_seed = -random_seed;
    }
#endif

}

double Diffusion::get_new_pos(const Walker* walker, int i, int j) {
    (void) walker;
    (void) i;
    (void) j;

#ifdef RNG_ZIG
    return DRanNormalZig32() * std;
#endif
#ifdef RNG_NUMREC
    return gaussian_deviate(&random_seed) * std;
#endif
}

double Diffusion::call_RNG() {
#ifdef RNG_ZIG
    return DRan_MWC8222();
#endif
#ifdef RNG_NUMREC
    return ran2(&random_seed);
#endif
}

void Diffusion::set_dt(double dt) {
    this->timestep = dt;
    this->std = sqrt(2 * D * dt);
}
