/* 
 * File:   DiTransform.cpp
 * Author: jorgmeister
 * 
 * Created on May 10, 2013, 1:46 PM
 */

#include "../../QMCheaders.h"

DiTransform::DiTransform(GeneralParams& gP, VariationalParams& vP)
: Orbitals(gP.n_p, gP.dim) {

    this->R = &(gP.R);

    gP.n_p /= 2;

    if (gP.system == "Diatom") {
        nucleus1 = new hydrogenicOrbitals(gP, vP);
        nucleus2 = new hydrogenicOrbitals(gP, vP);
    } else if (gP.system == "DoubleWell") {
        nucleus1 = new AlphaHarmonicOscillator(gP, vP);
        nucleus2 = new AlphaHarmonicOscillator(gP, vP);
    } else {
        std::cout << "System " << gP.system << " has not been implemented for DiTransform." << std::endl;
        exit(1);
    }

    gP.n_p *= 2;

    walker_nucleus1 = new Walker(n_p, dim);
    walker_nucleus2 = new Walker(n_p, dim);

}

void DiTransform::set_qnum_indie_terms(Walker* walker, int i) {

    walker->calc_r_i(i);

    walker_nucleus1->r.row(i) = walker->r.row(i);
    walker_nucleus2->r.row(i) = walker->r.row(i);

    walker_nucleus1->r(i, 0) += (*R) / 2;
    walker_nucleus2->r(i, 0) -= (*R) / 2;

    double shared = walker->get_r_i2(i) + 0.25 * (*R)*(*R);
    double comm_spec = walker->r(i, 0)*(*R);

    walker_nucleus1->r2(i) = shared + comm_spec;
    walker_nucleus2->r2(i) = shared - comm_spec;

    nucleus1->set_qnum_indie_terms(walker_nucleus1, i);
    nucleus2->set_qnum_indie_terms(walker_nucleus2, i);

}

double DiTransform::get_dell_alpha_phi(Walker* walker, int p, int q_num) {

    (void) walker;
    int sign = minusPower(q_num);

    return nucleus1->get_dell_alpha_phi(walker_nucleus1, p, q_num / 2)
            + sign * nucleus2->get_dell_alpha_phi(walker_nucleus2, p, q_num / 2);


}

double DiTransform::phi(const Walker* walker, int particle, int q_num) {

    (void) walker;
    int sign = minusPower(q_num);

    return nucleus1->phi(walker_nucleus1, particle, q_num / 2) +
            sign * nucleus2->phi(walker_nucleus2, particle, q_num / 2);
}

double DiTransform::del_phi(const Walker* walker, int particle, int q_num, int d) {

    (void) walker;
    int sign = minusPower(q_num);

    return nucleus1->del_phi(walker_nucleus1, particle, q_num / 2, d) +
            sign * nucleus2->del_phi(walker_nucleus2, particle, q_num / 2, d);
}

double DiTransform::lapl_phi(const Walker* walker, int particle, int q_num) {

    (void) walker;
    int sign = minusPower(q_num);

    return nucleus1->lapl_phi(walker_nucleus1, particle, q_num / 2) +
            sign * nucleus2->lapl_phi(walker_nucleus2, particle, q_num / 2);
}
