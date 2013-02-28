/* 
 * File:   Fokker_Planck.cpp
 * Author: jorgmeister
 * 
 * Created on October 12, 2012, 2:37 PM
 */

#include "../../QMCheaders.h"


Fokker_Planck::Fokker_Planck(int n_p, int dim, double timestep, seed_type random_seed, double D)
: Diffusion(n_p, dim, timestep, random_seed, D) {

}

double Fokker_Planck::get_g_ratio(const Walker* walker_post, const Walker* walker_pre) const {

    double g_ratio = 0;
    for (int i = 0; i < n_p; i++) {
        for (int j = 0; j < dim; j++) {
            g_ratio += (walker_pre->qforce(i,j) + walker_post->qforce(i,j))*
                    (D * timestep * 0.5 * (walker_pre->qforce(i,j) - walker_post->qforce(i,j))
                    - walker_post->r(i,j) + walker_pre->r(i,j));
        }
    }

    return exp(0.5*g_ratio);
}


