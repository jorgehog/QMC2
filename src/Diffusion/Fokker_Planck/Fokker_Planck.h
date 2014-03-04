#pragma once

#include "../Diffusion.h"
#include "../../Walker/Walker.h"

/*! \brief Anisotropic diffusion by the Fokker-Planck equation.
 */
class Fokker_Planck : public Diffusion {
public:
    Fokker_Planck(int n_p, int dim, double timestep, seed_type random_seed, double D = 0.5);

    double get_g_ratio(const Walker* walker_post, const Walker* walker_pre) const;

    double get_new_pos(const Walker* walker, int i, int j) {
        return D * timestep * walker->qforce(i, j) + Diffusion::get_new_pos(walker, i, j);
    }

};
