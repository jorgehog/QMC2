/* 
 * File:   DiAtomic.cpp
 * Author: jorgmeister
 * 
 * Created on May 10, 2013, 1:46 PM
 */

#include "../../QMCheaders.h"

DiAtomic::DiAtomic(GeneralParams& gP, VariationalParams& vP)
: Orbitals(gP.n_p, gP.dim) {

    this->R = &(gP.R);

//    vP.alpha = getAlpha(vP, n_p / 2.0);

    nucleus1 = new hydrogenicOrbitals(gP, vP, gP.n_p / 2);
    nucleus2 = new hydrogenicOrbitals(gP, vP, gP.n_p / 2);

    walker_nucleus1 = new Walker(n_p, dim);
    walker_nucleus2 = new Walker(n_p, dim);

}

double DiAtomic::getAlpha(VariationalParams& vP, double k0) {

    double tolerance = 1E-5;
    int NMAX = 10000;

    double k_pre;

    double k_new = k0;
    int N = 0;

    while (fabs(k_pre - k_new) > tolerance && N < NMAX) {
        k_pre = k_new;
        k_new = 1 + exp(-(*R) / k_pre);

        N++;
    }

    std::cout << "Found optimal alpha at " << k_new << "  " << 2 * k_new / n_p << std::endl;

    return k_new;
    //    return 2 * k_new / n_p;

}

void DiAtomic::set_qnum_indie_terms(Walker* walker, int i) {

    walker->calc_r_i(i);

    walker_nucleus1->r(i, 0) = walker_nucleus2->r(i, 0) = walker->r(i, 0);
    walker_nucleus1->r(i, 1) = walker_nucleus2->r(i, 1) = walker->r(i, 1);
    walker_nucleus1->r(i, 2) = walker_nucleus2->r(i, 2) = walker->r(i, 2);

    walker_nucleus1->r(i, 0) += (*R) / 2;
    walker_nucleus2->r(i, 0) -= (*R) / 2;

    double shared = walker->get_r_i2(i) + 0.25 * (*R)*(*R);
    double comm_spec = walker->r(i, 0)*(*R);

    walker_nucleus1->r2(i) = shared + comm_spec;
    walker_nucleus2->r2(i) = shared - comm_spec;

    nucleus1->set_qnum_indie_terms(walker_nucleus1, i);
    nucleus2->set_qnum_indie_terms(walker_nucleus2, i);

}

double DiAtomic::get_variational_derivative(const Walker* walker, int n) {
    
    walker_nucleus1->phi = walker_nucleus2->phi = walker->phi;
    walker_nucleus1->inv = walker_nucleus2->inv = walker->inv;
        
    return nucleus1->get_variational_derivative(walker_nucleus1, n) +
            nucleus2->get_variational_derivative(walker_nucleus2, n);
}

double DiAtomic::phi(const Walker* walker, int particle, int q_num) {
   
    (void) walker;
    
    return nucleus1->phi(walker_nucleus1, particle, q_num) +
            nucleus2->phi(walker_nucleus2, particle, q_num);
}

double DiAtomic::del_phi(const Walker* walker, int particle, int q_num, int d) {
    return nucleus1->del_phi(walker_nucleus1, particle, q_num, d) +
            nucleus2->del_phi(walker_nucleus2, particle, q_num, d);
}

double DiAtomic::lapl_phi(const Walker* walker, int particle, int q_num) {
    return nucleus1->lapl_phi(walker_nucleus1, particle, q_num) +
            nucleus2->lapl_phi(walker_nucleus2, particle, q_num);
}
