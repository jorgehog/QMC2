/* 
 * File:   HartreeFock.cpp
 * Author: jorgmeister
 * 
 * Created on April 4, 2013, 2:46 PM
 */

#include <armadillo>

class QMC;
class BasisFunctions;
class Walker;

#include "../Orbitals/Orbitals.h"
#include "HartreeFock.h"
#include <iostream>
#include <assert.h>

HartreeFock::HartreeFock(int m, Orbitals* sp_basis, double thresh) {

    this->thresh = 1E-6;

    assert(m / 2 <= sp_basis->max_implemented);

    this->n_p = sp_basis->n_p; //friends <3
    this->m = m;
    this->m2 = m / 2;

    this->sp_basis = sp_basis;

    H_hf.zeros(m, m);

    C.eye(m, m);
    e_HF.zeros(m);

    e_sp.zeros(m2);
    V.set_size(m2, m2);

    for (int i = 0; i < m2; i++) {
        for (int j = 0; j < m2; j++) {
            V(i, j) = arma::mat(m2, m2);
        }
    }

}

void HartreeFock::run_method() {

    setup_sp_e();
    setup_V();

    double e_new = 0;
    double e_old;

    std::cout << e_new << std::endl;
    do {
        e_old = e_new;

        make_Hamiltonian();

        std::cout << H_hf << std::endl;

        arma::eig_sym(e_HF, C, H_hf);
        C = C.st();
        sort();

        e_new = e_HF(0);

        std::cout << e_HF.st() << std::endl;

    } while (fabs(e_old - e_new) > thresh);

    double E_hf = calc_HF_E();

    std::cout << "fin " << E_hf << std::endl;
    C.resize(n_p, m);
    std::cout << C << std::endl;

    double E = 0;
    for (int p = 0; p < n_p; p++) {
        for (int i = 0; i < m; i++) {
                E += C(p, i) * C(p, i) * e_sp(i / 2);
        }
    }
    
    for (int a = 0; a < n_p; a++) {
        for (int b = 0; b < n_p; b++) {
            for (int c = 0; c < n_p; c++) {
                for (int d = 0; d < n_p; d++) {
                    for (int i = 0; i < m; i++) {
                        for (int j = 0; j < m; j++) {
                            for (int k = 0; k < m; k++) {
                                for (int l = 0; l < m; l++) {
                                    E += 0.25*C(a, i)*C(b, j)*C(c, k)*C(d, l)*get_V(i, j, k, l);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    std::cout << "energy: " << E << std::endl;
    V.reset();
    H_hf.reset();
    e_sp.reset();



}

void HartreeFock::sort() {

    using namespace arma;

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {

            if (C(i, j) < 0) {
                C(i, j) = -C(i, j);
            }

            if (C(i, j) < 1E-2) {
                C(i, j) = 0;
            }
        }
    }

    uvec indices = stable_sort_index(e_HF);

    vec e_HF_unsorted = e_HF;
    mat C_unsorted = C;


    for (int i = 0; i < m; i++) {
        C.col(i) = C_unsorted.col(indices(i));
        e_HF(i) = e_HF_unsorted(indices(i));
    }

    e_HF_unsorted.reset();
    C_unsorted.reset();

}

double HartreeFock::calc_HF_E() {
    int a, b, g, d, i, j;

    double E_hf = 0;

    for (i = 0; i < n_p; i++) {

        for (a = 0; a < m; a++) {
            E_hf += C(i, a) * C(i, a) * e_sp(a / 2);

            for (j = 0; j < n_p; j++) {

                for (b = 0; b < m; b++) {
                    for (g = 0; g < m; g++) {
                        for (d = 0; d < m; d++) {
                            E_hf += C(i, a) * C(j, b) * C(i, g) * C(j, d) * get_V(a, b, g, d);
                        }
                    }
                }
            }
        }
    }
}

void HartreeFock::make_Hamiltonian() {

    int a, g, b, d, i;
    using namespace std;
    H_hf.zeros();
    for (a = 0; a < m; a++) {

        H_hf(a, a) = e_sp(a / 2);

        for (g = 0; g < m; g++) {
            for (i = 0; i < n_p; i++) {
                for (b = 0; b < m; b++) {
                    for (d = 0; d < m; d++) {
                        H_hf(a, g) += C(i, b) * C(i, d) * get_V(a, b, g, d);
                    }
                }
            }
        }
    }
}

double HartreeFock::get_V(int a, int b, int g, int d) {

    int s0 = a % 2;
    int s1 = b % 2;
    int s2 = g % 2;
    int s3 = d % 2;

    bool s_conserved_dir = (s0 == s2) && (s1 == s3);
    bool s_conserved_exc = (s0 == s3) && (s1 == s2);

    double v1 = 0;
    double v2 = 0;

    if (s_conserved_dir) {
        v1 = V(a / 2, b / 2)(g / 2, d / 2);
    }

    if (s_conserved_exc) {
        v2 = V(a / 2, b / 2)(d / 2, g / 2);
    }

    //    std::cout << "--" << a << "--" << b << "--" << g << "--" << d << "--" << std::endl;
    //    std::cout << v1 << "   " << v2 << std::endl;
    return v1 - v2;
}

void HartreeFock::setup_V() {

    int a, b, g, d;
    arma::uvec qnum_set(4);

    for (a = 0; a < m2; a++) {
        qnum_set(0) = a;

        for (b = 0; b < m2; b++) {
            qnum_set(1) = b;

            for (g = a; g < m2; g++) {
                qnum_set(2) = g;

                for (d = b; d < m2; d++) {
                    qnum_set(3) = d;

                    V(a, b)(g, d) =
                            V(g, b)(a, d) =
                            V(a, d)(g, b) =
                            V(g, d)(a, b) =
                            sp_basis->get_coulomb_element(qnum_set);

                }
            }
        }
    }
}

void HartreeFock::setup_sp_e() {

    for (int a = 0; a < m2; a++) {
        e_sp(a) = sp_basis->get_sp_energy(a);
    }
}