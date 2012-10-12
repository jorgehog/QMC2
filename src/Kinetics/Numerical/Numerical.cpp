/* 
 * File:   Numerical.cpp
 * Author: jorgmeister
 * 
 * Created on October 12, 2012, 2:39 PM
 */

#include "../../QMCheaders.h"


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

