#pragma once

#include "../Diffusion.h"


namespace QMC2
{


/*! \brief Simple Isotropic diffusion model.
 */
class Simple : public Diffusion {
public:
    Simple(int n_p, int dim, double timestep, seed_type random_seed, double D = 0.5);

    double get_new_pos(const Walker* walker, int i, int j) {
        return Diffusion::get_new_pos(walker, i, j);
    }

    double get_g_ratio(const Walker* walker_post, const Walker* walker_pre) const {
        (void) walker_post;
        (void) walker_pre;

        return 1;
    }

};

}
