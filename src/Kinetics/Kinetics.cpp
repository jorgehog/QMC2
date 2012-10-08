/* 
 * File:   Kinetics.cpp
 * Author: jorgehog
 * 
 * Created on 13. april 2012, 17:45
 */

#include "../QMCheaders.h"

using namespace std;

Kinetics::Kinetics() {

}

Kinetics::Kinetics(int n_p, int dim) {
    this->n_p = n_p;
    this->n2 = n_p / 2;
    this->dim = dim;
}

NoKinetics::NoKinetics() {

}

Numerical::Numerical() {

}

Numerical::Numerical(GeneralParams & gP)
: Kinetics(gP.n_p, gP.dim) {
    this->h = gP.h;
    this->h2 = 1. / (h * h);

    wfminus = new Walker(n_p, dim);
    wfplus = new Walker(n_p, dim);

}

double Numerical::get_KE(const Walker* walker) {

    double e_kinetic, wf, wf_min, wf_plus;

    wf = walker->value;

    //    qmc->get_wf_value(walker); FUNKER
    //    if (abs(wf - walker->value) > 0.001) cout << "fail" << endl;

    for (int i = 0; i < n_p; i++) {
        for (int j = 0; j < dim; j++) {
            wfplus->r(i, j) = wfminus->r(i, j) = walker->r(i, j);
        }
    }

    e_kinetic = 0;
    for (int i = 0; i < n_p; i++) {
        for (int j = 0; j < dim; j++) {
            wfplus->r(i, j) = walker->r(i, j) + h;
            wfminus->r(i, j) = walker->r(i, j) - h;

            wfplus->make_rel_matrix();
            wfminus->make_rel_matrix();

            wfplus->calc_r_i2();
            wfminus->calc_r_i2();

            qmc->get_wf_value(wfplus);
            qmc->get_wf_value(wfminus);

            wf_min = wfminus->value;
            wf_plus = wfplus->value;

            e_kinetic -= (wf_min + wf_plus - 2 * wf);
            wfplus->r(i, j) = walker->r(i, j);
            wfminus->r(i, j) = walker->r(i, j);
        }
    }

    e_kinetic = 0.5 * h2 * e_kinetic / wf;

    return e_kinetic;
}

void Numerical::get_QF(Walker* walker) {
    int i, j;
    double wf_min, wf_plus, wf;

    wf = walker->value;

    for (i = 0; i < n_p; i++) {
        for (j = 0; j < dim; j++) {
            wfplus->r(i, j) = wfminus->r(i, j) = walker->r(i, j);
        }
    }


    for (i = 0; i < n_p; i++) {
        for (j = 0; j < dim; j++) {
            wfplus->r(i, j) = walker->r(i, j) + h;
            wfminus->r(i, j) = walker->r(i, j) - h;

            wfplus->make_rel_matrix();
            wfminus->make_rel_matrix();

            wfplus->calc_r_i2();
            wfminus->calc_r_i2();

            qmc->get_wf_value(wfplus);
            qmc->get_wf_value(wfminus);

            wf_min = wfminus->value;
            wf_plus = wfplus->value;

            walker->qforce(i, j) = (wf_plus - wf_min) / (wf * h);

            wfplus->r(i, j) = walker->r(i, j);
            wfminus->r(i, j) = walker->r(i, j);

        }
    }
}

void Numerical::get_necessities_IS(Walker* walker) const {
    qmc->get_wf_value(walker);
}

void Numerical::update_necessities_IS(const Walker* walker_pre, Walker* walker_post, int particle) const {
    qmc->get_wf_value(walker_post);
}

void Numerical::calculate_energy_necessities(Walker* walker) const {
    //No necessities.
}

double Numerical::get_spatial_ratio_IS(const Walker* walker_post, const Walker* walker_pre, int particle) const {
    return walker_post->value / walker_pre->value;
}

void Numerical::update_walker_IS(Walker* walker_pre, const Walker* walker_post, int particle) const {
    walker_pre->value = walker_post->value;
}

void Numerical::copy_walker_BF(const Walker* parent, Walker* child) const {
    child->value = parent->value;
}

void Numerical::copy_walker_IS(const Walker* parent, Walker* child) const {
    //nothing to copy. Might need for DMC IS
}

void Numerical::reset_walker_IS(const Walker* walker_pre, Walker* walker_post, int particle) const {
    //Nothing to reset;
}

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
            walker_pre->spatial_grad(i, j) = walker_post->spatial_grad(i, j);
        }
    }

    for (int i = 0; i < n_p; i++) {
        for (int j = 0; j < dim; j++) {
            walker_pre->jast_grad(i, j) = walker_post->jast_grad(i, j);
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
            xterm += walker->jast_grad(i, j) * walker->spatial_grad(i, j);
        }
    }

    e_kinetic = 2 * xterm + walker->lapl_sum;


    return -0.5 * e_kinetic;
}

void Closed_form::get_QF(Walker* walker) {
    int i, j;

    for (i = 0; i < n_p; i++) {
        for (j = 0; j < dim; j++) {
            walker->qforce(i, j) = 2 * (walker->jast_grad(i, j) + walker->spatial_grad(i, j));
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
    qmc->get_gradients(walker_post, particle);
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
            walker_post->spatial_grad(i, j) = walker_pre->spatial_grad(i, j);
        }
    }
}
