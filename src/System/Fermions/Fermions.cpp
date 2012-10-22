/* 
 * File:   Fermions.cpp
 * Author: jorgmeister
 * 
 * Created on October 12, 2012, 2:44 PM
 */

#include "../../QMCheaders.h"

Fermions::Fermions(GeneralParams & gP, Orbitals* orbital)
: System(gP.n_p, gP.dim, orbital) {

    I = arma::zeros<arma::rowvec > (n_p);

}

void Fermions::get_spatial_grad_full(Walker* walker) const {
    using namespace arma;

    for (int i = 0; i < n_p; i++) {
        walker->spatial_grad(i, span()) = strans(walker->dell_phi(i) * walker->inv(span(), i));
    }

}

void Fermions::get_spatial_grad(Walker* walker, int particle) const {
    using namespace arma;

    for (int i = start; i < n2 + start; i++) {
        walker->spatial_grad(i, span()) = strans(walker->dell_phi(i) * walker->inv(span(), i));
    }

}

void Fermions::make_merged_inv(Walker* walker) {
    using namespace arma;

    int _end;
    for (int _start = 0; _start < n_p; _start += n2) {
        _end = n2 + _start - 1;
        walker->inv(span(0, n2 - 1), span(_start, _end)) = arma::inv(walker->phi(span(_start, _end), span(0, n2 - 1)));
    }
}

double Fermions::get_spatial_ratio(const Walker* walker_pre, const Walker* walker_post, int particle) const {
    int q_num;
    double s_ratio;

    s_ratio = 0;
    for (q_num = 0; q_num < n2; q_num++) {
        s_ratio += walker_post->phi(particle, q_num) * walker_pre->inv(q_num, particle);
    }

    return s_ratio;
}

void Fermions::update_inverse(const Walker* walker_old, Walker* walker_new, int particle) {
    int k, l, j;
    double sum;


    for (j = start; j < n2 + start; j++) {
        if (j == particle) {
            I(j) = walker_new->spatial_ratio;
        } else {
            sum = 0;
            for (l = 0; l < n2; l++) {
                sum += walker_new->phi(particle, l) * walker_old->inv(l, j);
            }
            I(j) = sum;
        }
    }

    for (k = 0; k < n2; k++) {
        for (j = start; j < n2 + start; j++) {
            if (j == particle) {
                walker_new->inv(k, j) = walker_old->inv(k, particle) / walker_new->spatial_ratio;
            } else {
                walker_new->inv(k, j) = walker_old->inv(k, j) - walker_old->inv(k, particle) / walker_new->spatial_ratio * I(j);
            }
        }
    }
}

double Fermions::get_spatial_lapl_sum(const Walker* walker) const {
    int i, j;
    double sum;

    sum = 0;
    for (i = 0; i < n_p; i++) {
        for (j = 0; j < n_p / 2; j++) {
            sum += orbital->lapl_phi(walker, i, j) * walker->inv(j, i);
        }
    }

    return sum;
}
