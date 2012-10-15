/* 
 * File:   Importance.cpp
 * Author: jorgmeister
 * 
 * Created on October 12, 2012, 2:43 PM
 */

#include "../../QMCheaders.h"


Importance::Importance(GeneralParams & gP)
: Sampling(gP.n_p, gP.dim) {
    diffusion = new Fokker_Planck(n_p, dim, 0, gP.random_seed, gP.D);
    
}

void Importance::calculate_energy_necessities_CF(Walker* walker) const {
    //CF IS has no spesific necessities
}

void Importance::copy_walker(const Walker* parent, Walker* child) const {
    qmc->get_kinetics_ptr()->copy_walker_IS(parent, child);

    child->qforce = parent->qforce;
}

void Importance::update_necessities(const Walker* walker_pre, Walker* walker_post, int particle) {
    qmc->get_kinetics_ptr()->update_necessities_IS(walker_pre, walker_post, particle);
    qmc->get_kinetics_ptr()->get_QF(walker_post);
}

void Importance::update_walker(Walker* walker_pre, const Walker* walker_post, int particle) const {

    qmc->get_kinetics_ptr()->update_walker_IS(walker_pre, walker_post, particle);

    for (int i = 0; i < n_p; i++) {
        for (int j = 0; j < dim; j++) {
            walker_pre->qforce(i,j) = walker_post->qforce(i,j);
        }
    }
}

void Importance::reset_walker(const Walker* walker_pre, Walker* walker_post, int particle) const {
    qmc->get_kinetics_ptr()->reset_walker_IS(walker_pre, walker_post, particle);
}

double Importance::get_spatial_ratio(const Walker* walker_post, const Walker* walker_pre, int particle) const {
    return qmc->get_kinetics_ptr()->get_spatial_ratio_IS(walker_post, walker_pre, particle);
}

void Importance::get_necessities(Walker* walker) {
    qmc->get_kinetics_ptr()->get_necessities_IS(walker);
    qmc->get_kinetics_ptr()->get_QF(walker);
}
