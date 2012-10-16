/* 
 * File:   Fermions.cpp
 * Author: jorgmeister
 * 
 * Created on October 12, 2012, 2:44 PM
 */

#include "../../QMCheaders.h"

Fermions::Fermions(GeneralParams & gP, Orbitals* orbital)
: System(gP.n_p, gP.dim, orbital) {

    n2 = n_p / 2;

    s_up = arma::zeros<arma::mat > (n_p / 2, n_p / 2);
    s_down = arma::zeros<arma::mat > (n_p / 2, n_p / 2);
}

void Fermions::initialize(Walker* walker) {
    make_merged_inv(walker);
}

void Fermions::get_spatial_grad(Walker* walker, int particle) const {
    int i, j, k, start;


    start = n2 * (particle >= n2);

//    walker->spatial_grad(span(start, start + n2 - 1), span(0, dim - 1)) = arma::zeros<mat > (n2, dim);

    //updating the part with the same spin as the moved particle
//    cout << walker->dell_phi << endl;
//    cout << "-----" << endl;
    for (i = start; i < n2 + start; i++) {
        walker->spatial_grad(i, span()) = strans(walker->dell_phi(i)*walker->inv(span(), i));
    }

}

void Fermions::initialize_slaters(const Walker* walker) {

    for (int i = 0; i < n2; i++) {
        for (int q_num = 0; q_num < n2; q_num++) {
            s_up(i, q_num) = orbital->phi(walker, i, q_num);
            s_down(i, q_num) = orbital->phi(walker, i + n2, q_num);
        }
    }

}

double Fermions::get_det() {
    return arma::det(s_up) * arma::det(s_down);
}

void Fermions::invert_slaters() {
    s_up = arma::inv(s_up);
    s_down = arma::inv(s_down);

}

void Fermions::make_merged_inv(Walker* walker) {
//    int i, j;

    //    initialize_slaters(walker);
    //    invert_slaters();
    //
    //    //merging the inverse matrices
    //    for (i = 0; i < n2; i++) {
    //        for (j = 0; j < n2; j++) {
    //            walker->inv(i, j) = s_up(i, j);
    //            walker->inv(i, j + n2) = s_down(i, j);
    //        }
    //    }
    int end;
    for (int start = 0; start < n_p; start += n2) {
        end = n2 + start - 1;
        walker->inv(span(0, n2 - 1), span(start, end)) = arma::inv(walker->phi(span(start, end), span(0, n2 - 1)));
    }
}

double Fermions::get_spatial_wf(const Walker* walker) {
    initialize_slaters(walker);
    return get_det();
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

void Fermions::update_inverse(const Walker* walker_old, Walker* walker_new, int particle) const {
    int k, l, j, start;
    double sum;

    start = n2 * (particle >= n2);
    arma::rowvec I = arma::zeros<arma::rowvec > (n2);

    for (j = start; j < n2 + start; j++) {
        sum = 0;
        for (l = 0; l < n2; l++) {
            sum += walker_new->phi(particle, l) * walker_old->inv(l, j);
        }
        I(j - start) = sum;
    }

    //updating the part of the inverse with the same spin as particle i
    for (k = 0; k < n2; k++) {
        for (j = start; j < n2 + start; j++) {
            if (j == particle) {
                walker_new->inv(k, j) = walker_old->inv(k, particle) / walker_new->spatial_ratio;
            } else {

                walker_new->inv(k, j) = walker_old->inv(k, j) - walker_old->inv(k, particle) / walker_new->spatial_ratio * I(j - start);
            }
        }
    }
}

void Fermions::calc_for_newpos(const Walker* walker_old, Walker* walker_new, int particle) const {
    update_inverse(walker_old, walker_new, particle);
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

void Fermions::update_walker(Walker* walker_pre, const Walker* walker_post, int particle) const {
    int start = n2 * (particle >= n2);
    int end = start + n2 - 1;

    //Reseting the parts with the same spin as the moved particle
    walker_pre->inv(span(0, n2 - 1), span(start, end)) = walker_post->inv(span(0, n2 - 1), span(start, end));

}

void Fermions::copy_walker(const Walker* parent, Walker* child) const {
    child->inv = parent->inv;
}

void Fermions::reset_walker_ISCF(const Walker* walker_pre, Walker* walker_post, int particle) const {
    int start = n2 * (particle >= n2);
    int end = start + n2 - 1;

    //Reseting the part of the inverse with the same spin as particle i
    walker_post->inv(span(0, n2 - 1), span(start, end)) = walker_pre->inv(span(0, n2 - 1), span(start, end));

}