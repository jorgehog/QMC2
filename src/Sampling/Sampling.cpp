/* 
 * File:   Sampling.cpp
 * Author: jorgehog
 * 
 * Created on 15. juni 2012, 18:44
 */

#include "../QMCheaders.h"

Sampling::Sampling(int n_p, int dim) {
    this->n_p = n_p;
    this->n2 = n_p / 2;
    this->dim = dim;
}

Sampling::Sampling() {

}

void Sampling::set_trial_pos(Walker* walker, bool set_pos) {

    if (set_pos) {
        for (int i = 0; i < n_p; i++) {
            for (int j = 0; j < dim; j++) {
                walker->r(i, j) = diffusion->Diffusion::get_new_pos(walker, i, j);
            }
        }
    }
    
    walker->calc_r_i2();

    set_trial_states(walker);

    walker->make_rel_matrix();
    
    qmc->get_system_ptr()->initialize(walker);

    qmc->get_jastrow_ptr()->get_dJ_matrix(walker);
    
    get_necessities(walker);

}

void Sampling::set_trial_states(Walker * walker) {

    for (int i = 0; i < n_p; i++) {
       
        qmc->get_orbitals_ptr()->set_qnum_indie_terms(walker, i);
     
        for (int j = 0; j < n2; j++) {
            walker->phi(i, j) = qmc->get_orbitals_ptr()->phi(walker, i, j);
        }
        for (int j = 0; j < n2; j++) {
            for (int k = 0; k < dim; k++) {
                walker->dell_phi(i)(k, j) = qmc->get_orbitals_ptr()->del_phi(walker, i, j, k);
            }
        }
    }
}

void Sampling::update_pos(const Walker* walker_pre, Walker* walker_post, int particle) const {

    for (int j = 0; j < dim; j++) {
        walker_post->r(particle, j) = walker_pre->r(particle, j)
                + diffusion->get_new_pos(walker_pre, particle, j);
    }

    for (int j = 0; j < n_p; j++) {
        if (j != particle) {
            walker_post->r_rel(particle, j) = walker_post->r_rel(j, particle)
                    = walker_post->calc_r_rel(particle, j);
        }
    }

    walker_post->calc_r_i2(particle);
    qmc->get_orbitals_ptr()->set_qnum_indie_terms(walker_post, particle);

    for (int j = 0; j < n2; j++) {
        walker_post->phi(particle, j) = qmc->get_orbitals_ptr()->phi(walker_post, particle, j);
    }

    for (int j = 0; j < n2; j++) {
        for (int k = 0; k < dim; k++) {
            walker_post->dell_phi(particle)(k, j) = qmc->get_orbitals_ptr()->del_phi(walker_post, particle, j, k);
        }
    }

    qmc->get_jastrow_ptr()->get_dJ_matrix(walker_post, particle);
    walker_post->spatial_ratio = qmc->get_system_ptr()->get_spatial_ratio(walker_pre, walker_post, particle);
    qmc->get_system_ptr()->calc_for_newpos(walker_pre, walker_post, particle);

    update_necessities(walker_pre, walker_post, particle);

}
