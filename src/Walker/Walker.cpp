/* 
 * File:   Walker.cpp
 * Author: jorgehog
 * 
 * Created on 30. mars 2012, 16:50
 */

#include "../QMCheaders.h"

using namespace arma;

Walker::Walker(int n_p, int dim, bool alive) {
    this->dim = dim;
    this->n_p = n_p;
    this->n2 = n_p / 2;


    if (alive) {
        is_murdered = false;
    } else {
        is_murdered = true;
    }

    r = zeros<mat > (n_p, dim);
    r_rel = zeros<mat > (n_p, n_p);
    qforce = zeros<mat > (n_p, dim);
    jast_grad = zeros<mat > (n_p, dim);
    spatial_grad = zeros<mat > (n_p, dim);

    r2 = zeros(1, n_p);

    value = 0;
    lapl_sum = 0;
    spatial_ratio = 0;
    inv = zeros<mat > (n2, n_p);

    phi = zeros<mat > (n_p, n2);
    dell_phi = arma::field<mat > (n_p, 1);
    for (int i = 0; i < n_p; i++){
        dell_phi(i) = zeros<mat>(dim, n2);
    }
    
    dJ = zeros<cube > (n_p, n_p, dim);

}

void Walker::calc_r_i2(int i) {

    double r2i = 0;

    for (int j = 0; j < dim; j++) {
        r2i += r(i, j) * r(i, j);
    }

    r2[i] = r2i;

}

void Walker::calc_r_i2() {

    for (int i = 0; i < n_p; i++) {
        this->calc_r_i2(i);
    }

}

double Walker::abs_relative(int i, int j) const {
    int k;
    double r_ij, tmp;

    r_ij = 0;
    for (k = 0; k < dim; k++) {
        tmp = (r(i, k) - r(j, k));
        r_ij += tmp*tmp;
    }
    r_ij = sqrt(r_ij);

    return r_ij;
}

double Walker::get_r_i2(int i) const {
    return r2(i);
}

void Walker::make_rel_matrix() {
    int i, j;


    for (i = 0; i < n_p - 1; i++) {
        for (j = i + 1; j < n_p; j++) {
            r_rel(i, j) = r_rel(j, i) = abs_relative(i, j);
        }
    }
}

bool Walker::is_singular() const {
    int i, j;
    double eps;

    eps = 0.1;

    for (i = 0; i < n_p; i++) {
        for (j = i + 1; j < n_p; j++) {
            if (r_rel(i, j) < eps) {
                return true;
            }
        }
    }

    return false;
}

bool Walker::check_bad_qforce() {
    double qforce_test = 0;

    for (int i = 0; i < n_p; i++) {
        for (int j = 0; j < dim; j++) {
            double tmp = qforce(i, j);
            if (tmp * tmp > qforce_test) {
                qforce_test = tmp*tmp;
            }
        }
    }

    //expression from regression table
    if (qforce_test > (500 * n_p - 6000)*(n_p >= 12) + 1000) {
        return true;
    } else {
        return false;
    }

}

void Walker::print(std::string header) {
    using namespace std;

    cout << endl;
    cout << "---- ---- ---- " << header << " ---- ---- ----" << endl;

    cout << "r\n" << this->r << endl;
    cout << "r_rel\n" << this->r_rel << endl;
    cout << "r2\n" << this->r2 << endl;

    cout << "S grad\n" << this->spatial_grad << endl;
    cout << "J grad\n" << this->jast_grad << endl;

    cout << "Qforce\n" << this->qforce << endl;
    cout << "inv\n" << this->inv << endl;

    cout << "Lapl_sum\t" << this->lapl_sum << endl;
    cout << "Value\t" << this->value << endl;
    cout << "S ratio\t" << this->spatial_ratio << endl;

    cout << "---- ---- ---- ---- ---- ---- ----" << endl;
    cout << endl;
}