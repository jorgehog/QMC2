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
    system->get_spatial_grad(walker, 0);
    system->get_spatial_grad(walker, n2);
}

void QMC::get_laplsum(Walker* walker) const {
    walker->lapl_sum = system->get_spatial_lapl_sum(walker) + jastrow->get_lapl_sum(walker);
}

void QMC::get_wf_value(Walker* walker) const {
    walker->value = system->get_spatial_wf(walker) * jastrow->get_val(walker);
}

void QMC::update_pos(const Walker* walker_pre, Walker* walker_post, int particle) const {

    for (int j = 0; j < dim; j++) {
        walker_post->r(particle, j) = walker_pre->r(particle, j)
                + sampling->get_new_pos(walker_pre, particle, j);
    }

    for (int j = 0; j < n_p; j++) {
        if (j != particle) {
            walker_post->r_rel(particle, j) = walker_post->r_rel(j, particle)
                    = walker_post->abs_relative(particle, j);
        }
    }

    walker_post->calc_r_i2(particle);
    for (int j = 0; j < n2; j++) {
        walker_post->phi(particle, j) = get_orbitals_ptr()->phi(walker_post, particle, j);

        for (int k = 0; k < dim; k++) {
            walker_post->dell_phi(particle)(k, j) = get_orbitals_ptr()->del_phi(walker_post, particle, j, k);
        }
    }

    jastrow->get_dJ_matrix(walker_post, particle);
    walker_post->spatial_ratio = system->get_spatial_ratio(walker_pre, walker_post, particle);
    system->calc_for_newpos(walker_pre, walker_post, particle);
    
    std::cout << "pre: "<<1-arma::accu(walker_pre->inv*walker_pre->phi)/n_p << std::endl;
    std::cout << "post: "<<1-arma::accu(walker_post->inv*walker_post->phi)/n_p << std::endl;
    sampling->update_necessities(walker_pre, walker_post, particle);

}

double QMC::get_acceptance_ratio(const Walker* walker_pre, const Walker* walker_post, int particle) const {
    double spatial_jast = sampling->get_spatialjast_ratio(walker_post, walker_pre, particle);
    double G = sampling->get_g_ratio(walker_post, walker_pre);

    return spatial_jast * spatial_jast * G;
}

void QMC::calculate_energy_necessities(Walker* walker) const {
    get_sampling_ptr()->calculate_energy_necessities_CF(walker);
    get_laplsum(walker);
    //    kinetics->calculate_energy_necessities(walker);
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

        update_pos(original, trial, particle);

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
    int i, j;

    for (i = 0; i < n_p; i++) {
        for (j = 0; j < dim; j++) {
            walker->qforce(i, j) = 2 * (walker->jast_grad(i, j) + walker->spatial_grad(i, j));
        }
    }
}