/* 
 * File:   QMC.cpp
 * Author: jorgehog
 * 
 * Created on 30. mars 2012, 17:42
 */

#include "../QMCheaders.h"

QMC::QMC(int n_p, int dim, int n_c,
        Sampling *sampling,
        System *system,
        Jastrow *jastrow) {

    this->n_p = n_p;
    this->dim = dim;
    this->n_c = n_c;
    this->n2 = n_p / 2;

    this->jastrow = jastrow;
    this->sampling = sampling;
    this->system = system;

    this->sampling->set_qmc_ptr(this);

    this->accepted = 0;

}

QMC::QMC() {

}

void QMC::add_output(OutputHandler* output_handler) {
    output_handler->set_qmc_ptr(this);
    this->output_handler.push_back(output_handler);
}

void QMC::dump_output() {

    for (std::vector<OutputHandler*>::iterator output_obj = output_handler.begin(); output_obj != output_handler.end(); ++output_obj) {
        (*output_obj)->dump();
    }

}

void QMC::finalize_output() {

    for (std::vector<OutputHandler*>::iterator output_obj = output_handler.begin(); output_obj != output_handler.end(); ++output_obj) {
        (*output_obj)->finalize();
    }

}

void QMC::get_gradients(const Walker* walker_pre, Walker* walker_post, int particle) const {
    jastrow->get_grad(walker_pre, walker_post, particle);
    system->get_spatial_grad(walker_post, particle);
}

void QMC::get_gradients(Walker* walker) const {
    jastrow->get_grad(walker);
    system->get_spatial_grad_full(walker);
}

double QMC::get_acceptance_ratio(const Walker* walker_pre, const Walker* walker_post, int particle) const {
    double spatial_jast = sampling->get_spatialjast_ratio(walker_post, walker_pre, particle);
    double G = sampling->get_g_ratio(walker_post, walker_pre);

    return spatial_jast * spatial_jast * G;
}

void QMC::set_spin_state(int particle) const {
    int start = n2 * (particle >= n2);
    int end = start + n2 - 1;
    system->set_spin_state(start, end);
    sampling->set_spin_state(start, end);
}

bool QMC::metropolis_test(double A) {
    double r = sampling->call_RNG();

    if (r <= A) {
        accepted++;
        return true;

    } else {
        return false;
    }
}

void QMC::calculate_energy_necessities(Walker* walker) const {
    sampling->calculate_energy_necessities(walker);
    get_laplsum(walker);
}

void QMC::update_walker(Walker* walker_pre, const Walker* walker_post, int particle) const {

    for (int i = 0; i < dim; i++) {
        walker_pre->r(particle, i) = walker_post->r(particle, i);
    }

    for (int i = 0; i < n_p; i++) {
        walker_pre->r_rel(i, particle) = walker_pre->r_rel(particle, i) = walker_post->r_rel(i, particle);

        for (int k = 0; k < dim; k++) {
            walker_pre->dJ(particle, i, k) = walker_post->dJ(particle, i, k);
            walker_pre->dJ(i, particle, k) = walker_post->dJ(i, particle, k);
        }
    }

    for (int i = 0; i < n2; i++) {
        walker_pre->phi(particle, i) = walker_post->phi(particle, i);
        walker_pre->dell_phi(particle)(arma::span(), i) = walker_post->dell_phi(particle)(arma::span(), i);
    }

    walker_pre->r2[particle] = walker_post->r2[particle];

    system->update_walker(walker_pre, walker_post, particle);
    sampling->update_walker(walker_pre, walker_post, particle);
}

void QMC::reset_walker(const Walker* walker_pre, Walker* walker_post, int particle) const {
    for (int i = 0; i < dim; i++) {
        walker_post->r(particle, i) = walker_pre->r(particle, i);
    }

    for (int i = 0; i < n_p; i++) {
        walker_post->r_rel(i, particle) = walker_post->r_rel(particle, i) = walker_pre->r_rel(i, particle);

        for (int k = 0; k < dim; k++) {
            walker_post->dJ(particle, i, k) = walker_pre->dJ(particle, i, k);
            walker_post->dJ(i, particle, k) = walker_pre->dJ(i, particle, k);
        }
    }

    for (int i = 0; i < n2; i++) {
        walker_post->phi(particle, i) = walker_pre->phi(particle, i);
        walker_post->dell_phi(particle)(arma::span(), i) = walker_pre->dell_phi(particle)(arma::span(), i);
    }

    walker_post->r2[particle] = walker_pre->r2[particle];

    system->reset_walker(walker_pre, walker_post, particle);
    sampling->reset_walker(walker_pre, walker_post, particle);
}

void QMC::diffuse_walker(Walker* original, Walker* trial) {
    for (int particle = 0; particle < n_p; particle++) {

        set_spin_state(particle);
        sampling->update_pos(original, trial, particle);

        double A = get_acceptance_ratio(original, trial, particle);

        if (move_autherized(A)) {
            update_walker(original, trial, particle);
        } else {
            reset_walker(original, trial, particle);
        }

    }
}

void QMC::copy_walker(const Walker* parent, Walker* child) const {

    child->r2 = parent->r2;
    child->r = parent->r;
    child->r_rel = parent->r_rel;

    child->phi = parent->phi;
    child->dell_phi = parent->dell_phi;
    child->dJ = parent->dJ;

    child->ressurect();
    child->set_E(parent->get_E());
    system->copy_walker(parent, child);
    sampling->copy_walker(parent, child);

}

double QMC::get_KE(const Walker* walker) const {
    int i, j;
    double xterm, e_kinetic;

    e_kinetic = xterm = 0;

    for (i = 0; i < n_p; i++) {
        for (j = 0; j < dim; j++) {
            xterm += walker->jast_grad(i, j) * walker->spatial_grad(i, j);
        }
    }

    e_kinetic = 2 * xterm + walker->lapl_sum;


    return -0.5 * e_kinetic;
}

void QMC::get_QF(Walker* walker) const { 
    walker->qforce = 2*(walker->jast_grad + walker->spatial_grad);
}