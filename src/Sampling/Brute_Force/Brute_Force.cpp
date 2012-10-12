/* 
 * File:   Brute_Force.cpp
 * Author: jorgmeister
 * 
 * Created on October 12, 2012, 2:43 PM
 */

#include "../../QMCheaders.h"


Brute_Force::Brute_Force(GeneralParams & gP)
: Sampling(gP.n_p, gP.dim) {
    diffusion = new Simple(n_p, dim, 0, gP.random_seed, gP.D);

}

void Brute_Force::update_walker(Walker* walker_pre, const Walker* walker_post, int particle) const {
    walker_pre->value = walker_post->value;
}

void Brute_Force::get_necessities(Walker* walker) {
    qmc->get_wf_value(walker);
}

void Brute_Force::update_necessities(const Walker* walker_pre, Walker* walker_post, int particle) {
    qmc->get_wf_value(walker_post);
}

void Brute_Force::reset_walker(const Walker* walker_pre, Walker* walker_post, int particle) const {
    //    walker_post->value = walker_pre->value; //will be overwritten.
}

void Brute_Force::calculate_energy_necessities_CF(Walker* walker) const {
    qmc->get_system_ptr()->initialize(walker);
    qmc->get_gradients(walker);
}

double Brute_Force::get_spatial_ratio(const Walker* walker_post, const Walker* walker_pre, int particle) const {
    return walker_post->value / walker_pre->value;
}

void Brute_Force::copy_walker(const Walker* parent, Walker* child) const {
    child->value = parent->value;
}
