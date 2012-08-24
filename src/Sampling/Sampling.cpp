/* 
 * File:   Sampling.cpp
 * Author: jorgehog
 * 
 * Created on 15. juni 2012, 18:44
 */

#include "../QMCheaders.h"

using namespace std;

Sampling::Sampling(int n_p, int dim) {
    this->n_p = n_p;
    this->dim = dim;
}

void Sampling::set_trial_pos(Walker* walker, bool load_VMC_dist, std::ifstream* file) const {
    int i, j;

    if (load_VMC_dist) {

        double pos;
        for (i = 0; i < n_p; i++) {
            for (j = 0; j < dim; j++) {
                *file >> pos;
                walker->r(i, j) = pos;
            }
        }
    } else {
        for (i = 0; i < n_p; i++) {
            for (j = 0; j < dim; j++) {
                walker->r(i, j) = diffusion->Diffusion::get_new_pos(walker, i, j);
            }
        }
    }

    for (i = 0; i < n_p; i++) {
        walker->calc_r_i2(i);
    }
    
    walker->make_rel_matrix();

    get_necessities(walker);

}

double Sampling::get_new_pos(const Walker* walker_pre, int particle, int j) const {
    return diffusion->get_new_pos(walker_pre, particle, j);
}

double Sampling::get_g_ratio(const Walker* walker_post, const Walker* walker_pre, int particle) const {
    return diffusion->get_g_ratio(walker_post, walker_pre, particle);
}

double Sampling::get_branching_Gfunc(Walker* walker_pre, Walker* walker_post, double E_T) const {
    return diffusion->get_GBfunc(walker_pre, walker_post, E_T);
}

Brute_Force::Brute_Force(int n_p, int dim, double timestep, long random_seed, double D) : Sampling(n_p, dim) {
    is_importance = false;
    diffusion = new Simple(n_p, dim, timestep, random_seed, D);
}

void Brute_Force::update_walker(Walker* walker_pre, const Walker* walker_post, int particle) const {
    walker_pre->value = walker_post->value;
}

void Brute_Force::get_necessities(Walker* walker) const {
    qmc->get_wf_value(walker);
}

void Brute_Force::update_necessities(const Walker* walker_pre, Walker* walker_post, int particle) const {
    qmc->get_wf_value(walker_post);
}

void Brute_Force::reset_walker(const Walker* walker_pre, Walker* walker_post, int particle) const {
    //Nothing to reset;
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

Importance::Importance(int n_p, int dim, double timestep, long random_seed, double D) : Sampling(n_p, dim) {
    is_importance = true;
    diffusion = new Fokker_Planck(n_p, dim, timestep, random_seed, D);
}

void Importance::calculate_energy_necessities_CF(Walker* walker) const {
    //CF IS has no spesific necessities
}

void Importance::copy_walker(const Walker* parent, Walker* child) const {
    qmc->get_kinetics_ptr()->copy_walker_IS(parent, child);

    child->qforce = parent->qforce;
}

void Importance::update_necessities(const Walker* walker_pre, Walker* walker_post, int particle) const {
    qmc->get_kinetics_ptr()->update_necessities_IS(walker_pre, walker_post, particle);
    qmc->get_kinetics_ptr()->get_QF(walker_post);
}

void Importance::update_walker(Walker* walker_pre, const Walker* walker_post, int particle) const {

    qmc->get_kinetics_ptr()->update_walker_IS(walker_pre, walker_post, particle);

    for (int i = 0; i < n_p; i++) {
        for (int j = 0; j < dim; j++) {
            walker_pre->qforce(i, j) = walker_post->qforce(i, j);
        }
    }
}

void Importance::reset_walker(const Walker* walker_pre, Walker* walker_post, int particle) const {
    qmc->get_kinetics_ptr()->reset_walker_IS(walker_pre, walker_post, particle);
}

double Importance::get_spatial_ratio(const Walker* walker_post, const Walker* walker_pre, int particle) const {
    return qmc->get_kinetics_ptr()->get_spatial_ratio_IS(walker_post, walker_pre, particle);
}

void Importance::get_necessities(Walker* walker) const {
    qmc->get_kinetics_ptr()->get_necessities_IS(walker);
    qmc->get_kinetics_ptr()->get_QF(walker);

    if (walker->check_bad_qforce()) {
        Sampling::set_trial_pos(walker);
    }
}


