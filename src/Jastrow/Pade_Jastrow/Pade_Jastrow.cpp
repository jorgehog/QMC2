/* 
 * File:   Pade_Jastrow.cpp
 * Author: jorgmeister
 * 
 * Created on October 12, 2012, 2:38 PM
 */

#include "../../QMCheaders.h"

Pade_Jastrow::Pade_Jastrow(GeneralParams & gP, VariationalParams & vP)
: Jastrow(gP.n_p, gP.dim) {

    active = true;

    set_parameter(vP.beta, 0);

    a = arma::zeros<arma::mat > (n_p, n_p);
}

void Pade_Jastrow::initialize() {
    int i, j;
    double a_sym, a_asym;

    if (dim == 2) {
        a_sym = 1. / 3;
        a_asym = 1.0;
    } else if (dim == 3) {
        a_sym = 1. / 4;
        a_asym = 1. / 2;
    } else {
        std::cout << "Unable to initialize Jastrow paremters: Unknown dimension" << std::endl;
        exit(1);
    }

    for (i = 0; i < n_p; i++) {
        for (j = 0; j < n_p; j++) {

            if ((j < n2 && i < n2) || (j >= n2 && i >= n2)) {
                a(i, j) = a_sym;
            } else {
                a(i, j) = a_asym;
            }
        }
    }
}

double Pade_Jastrow::get_val(const Walker* walker) const {
    int i, j;
    double arg;

    arg = 0;
    for (i = 0; i < n_p - 1; i++) {
        for (j = i + 1; j < n_p; j++) {
            arg += a(i, j) * walker->r_rel(i, j) / (1.0 + beta * walker->r_rel(i, j));
        }
    }

    return exp(arg);
}

void Pade_Jastrow::get_dJ_matrix(Walker* walker, int i) const {
    double b_ij, factor;

    for (int j = 0; j < i; j++) {
        b_ij = 1.0 + beta * walker->r_rel(i, j);
        factor = a(i, j) / (walker->r_rel(i, j) * b_ij * b_ij);
        for (int k = 0; k < dim; k++) {
            walker->dJ(i, j, k) = (walker->r(i, k) - walker->r(j, k)) * factor;
            walker->dJ(j, i, k) = -walker->dJ(i, j, k);
        }
    }

    for (int j = i + 1; j < n_p; j++) {
        b_ij = 1.0 + beta * walker->r_rel(i, j);
        factor = a(i, j) / (walker->r_rel(i, j) * b_ij * b_ij);
        for (int k = 0; k < dim; k++) {
            walker->dJ(i, j, k) = (walker->r(i, k) - walker->r(j, k)) * factor;
            walker->dJ(j, i, k) = -walker->dJ(i, j, k);
        }
    }
}

void Pade_Jastrow::get_grad(Walker* walker) const {
    double sum;

    for (int i = 0; i < n_p; i++) {
        for (int k = 0; k < dim; k++) {

            sum = 0;
            for (int j = 0; j < i; j++) {
                sum += walker->dJ(i, j, k);
            }

            for (int j = i + 1; j < n_p; j++) {
                sum += walker->dJ(i, j, k);
            }

            walker->jast_grad(i, k) = sum;
        }
    }

}

void Pade_Jastrow::get_grad(const Walker* walker_pre, Walker* walker_post, int p) const {
    double sum;

    for (int i = 0; i < p; i++) {
        for (int k = 0; k < dim; k++) {
            walker_post->jast_grad(i, k) = walker_pre->jast_grad(i, k) + walker_post->dJ(i, p, k) - walker_pre->dJ(i, p, k);
        }
    }

    for (int k = 0; k < dim; k++) {

        sum = 0;
        for (int j = 0; j < p; j++) {
            sum += walker_post->dJ(p, j, k);
        }

        for (int j = p + 1; j < n_p; j++) {
            sum += walker_post->dJ(p, j, k);
        }

        walker_post->jast_grad(p, k) = sum;
    }

    for (int i = p + 1; i < n_p; i++) {
        for (int k = 0; k < dim; k++) {
            walker_post->jast_grad(i, k) = walker_pre->jast_grad(i, k) + walker_post->dJ(i, p, k) - walker_pre->dJ(i, p, k);
        }
    }

}

double Pade_Jastrow::get_j_ratio(const Walker* walker_new, const Walker* walker_old, int i) const {
    int j;
    double j_ratio;

    j_ratio = 0;
    for (j = 0; j < n_p; j++) {
        j_ratio += a(i, j) * (walker_new->r_rel(i, j) / (1.0 + beta * walker_new->r_rel(i, j)) -
                walker_old->r_rel(i, j) / (1.0 + beta * walker_old->r_rel(i, j)));
    }


    return exp(j_ratio);
}

double Pade_Jastrow::get_lapl_sum(const Walker * walker) const {
    int k, j, d;
    double sum1, sum2, b_kj;

    sum1 = sum2 = 0;

    for (k = 0; k < n_p; k++) {
        for (j = k + 1; j < n_p; j++) {
            b_kj = 1 + beta * walker->r_rel(k, j);
            sum2 += a(k, j) * (1 - beta * walker->r_rel(k, j)) / (walker->r_rel(k, j) * b_kj * b_kj * b_kj);
        }

        for (d = 0; d < dim; d++) {
            sum1 += walker->jast_grad(k, d) * walker->jast_grad(k, d);
        }
    }

    return sum1 + 2 * sum2;
}

double Pade_Jastrow::get_variational_derivative(const Walker* walker, int n) {

    double dbeta = 0.0;

    for (int i = 0; i < n_p - 1; i++) {
        for (int j = i + 1; j < n_p; j++) {
            double arg = walker->r_rel(i, j) / (1.0 + beta * walker->r_rel(i, j));
            dbeta -= a(i, j) * arg * arg;
        }
    }

    return dbeta;

}

