/* 
 * File:   Sampling.cpp
 * Author: jorgehog
 * 
 * Created on 15. juni 2012, 18:44
 */

#include "../QMCheaders.h"

Sampling::Sampling(int n_p, int dim) {
    this->n_p = n_p;
    this->dim = dim;
}

Sampling::Sampling() {

}

void Sampling::set_trial_pos(Walker* walker, bool load_VMC_dist, std::ifstream* file) {
    int i, j;

    if (load_VMC_dist) {

        double pos;
        for (i = 0; i < n_p; i++) {
            for (j = 0; j < dim; j++) {
                *file >> pos;
                walker->r.at(i,j) = pos;
            }
        }
    } else {
        for (i = 0; i < n_p; i++) {
            for (j = 0; j < dim; j++) {
                walker->r.at(i,j) = diffusion->Diffusion::get_new_pos(walker, i, j);
            }
        }
    }

    walker->calc_r_i2();
    walker->make_rel_matrix();
    get_necessities(walker);

}

double Sampling::get_new_pos(const Walker* walker_pre, int particle, int j) const {
    return diffusion->get_new_pos(walker_pre, particle, j);
}

double Sampling::get_g_ratio(const Walker* walker_post, const Walker* walker_pre) const {
    return diffusion->get_g_ratio(walker_post, walker_pre);
}

double Sampling::get_branching_Gfunc(double E_x, double E_y, double E_T) const {
    return diffusion->get_GBfunc(E_x, E_y, E_T);
}
