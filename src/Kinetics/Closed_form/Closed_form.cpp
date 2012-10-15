/* 
 * File:   Closed_form.cpp
 * Author: jorgmeister
 * 
 * Created on October 12, 2012, 2:40 PM
 */

#include "../../QMCheaders.h"


Closed_form::Closed_form() {

}

Closed_form::Closed_form(GeneralParams & gP)
: Kinetics(gP.n_p, gP.dim) {

}

void Closed_form::update_walker_IS(Walker* walker_pre, const Walker* walker_post, int particle) const {
    int start = n2 * (particle >= n2);

    qmc->get_system_ptr()->update_walker(walker_pre, walker_post, particle);

    for (int i = start; i < start + n2; i++) {
        for (int j = 0; j < dim; j++) {
            walker_pre->spatial_grad(i,j) = walker_post->spatial_grad(i,j);
        }
    }

    for (int i = 0; i < n_p; i++) {
        for (int j = 0; j < dim; j++) {
            walker_pre->jast_grad(i,j) = walker_post->jast_grad(i,j);
        }
    }
}

double Closed_form::get_KE(const Walker* walker) {
    int i, j;
    double xterm, e_kinetic;


    e_kinetic = xterm = 0;

    //the X-term
    for (i = 0; i < n_p; i++) {
        for (j = 0; j < dim; j++) {
            xterm += walker->jast_grad(i,j) * walker->spatial_grad(i,j);
        }
    }

    e_kinetic = 2 * xterm + walker->lapl_sum;


    return -0.5 * e_kinetic;
}

void Closed_form::get_QF(Walker* walker) {
    int i, j;

    for (i = 0; i < n_p; i++) {
        for (j = 0; j < dim; j++) {
            walker->qforce(i,j) = 2 * (walker->jast_grad(i,j) + walker->spatial_grad(i,j));
        }
    }
}

double Closed_form::get_spatial_ratio_IS(const Walker* walker_post, const Walker* walker_pre, int particle) const {
    return walker_post->spatial_ratio * qmc->get_jastrow_ptr()->get_j_ratio(walker_post, walker_pre, particle);
}

void Closed_form::get_necessities_IS(Walker* walker) const {
    qmc->get_system_ptr()->initialize(walker);
    qmc->get_gradients(walker);
}

void Closed_form::calculate_energy_necessities(Walker* walker) const {
    qmc->get_sampling_ptr()->calculate_energy_necessities_CF(walker);
    qmc->get_laplsum(walker);
}

void Closed_form::update_necessities_IS(const Walker* walker_pre, Walker* walker_post, int particle) const {
    walker_post->spatial_ratio = qmc->get_system_ptr()->get_spatial_ratio(walker_pre, walker_post, particle);
    qmc->get_system_ptr()->calc_for_newpos(walker_pre, walker_post, particle);
    qmc->get_gradients(walker_pre, walker_post, particle);
}

void Closed_form::copy_walker_BF(const Walker* parent, Walker* child) const {
    child->value = parent->value;
}

void Closed_form::copy_walker_IS(const Walker* parent, Walker* child) const {
    child->jast_grad = parent->jast_grad;
    child->spatial_grad = parent->spatial_grad;

    qmc->get_system_ptr()->copy_walker(parent, child);
}

void Closed_form::reset_walker_IS(const Walker* walker_pre, Walker* walker_post, int particle) const {

    qmc->get_system_ptr()->reset_walker_ISCF(walker_pre, walker_post, particle);

    int start = n2 * (particle >= n2);

    //updating the part with the same spin as the moved particle
    for (int i = start; i < n2 + start; i++) {
        for (int j = 0; j < dim; j++) {
            walker_post->spatial_grad(i,j) = walker_pre->spatial_grad(i,j);
        }
    }
}
