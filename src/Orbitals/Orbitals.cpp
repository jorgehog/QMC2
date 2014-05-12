/* 
 * File:   Orbitals.cpp
 * Author: jorgehog
 * 
 * Created on 26. juni 2012, 17:41
 */

#include "Orbitals.h"

#include "../Walker/Walker.h"
#include "../QMC/QMC.h"
#include "../BasisFunctions/BasisFunctions.h"

using namespace QMC2;

Orbitals::Orbitals(int n_p, int dim) {
    this->n_p = n_p;
    this->n2 = ceil(n_p / 2.0);
    this->dim = dim;

    h = 1E-4;
    h2 = 1 / (h * h);
    two_h = 1 / (2 * h);


    max_implemented = 28; //for 56 particles 
    basis_functions = new BasisFunctions*[max_implemented];

    dell_basis_functions = new BasisFunctions**[dim];
    for (int i = 0; i < dim; i++) {
        dell_basis_functions[i] = new BasisFunctions*[max_implemented];
    }

    lapl_basis_functions = new BasisFunctions*[max_implemented];

    nCap = n_p;
}

Orbitals::Orbitals() {

}

double Orbitals::num_diff(const Walker* walker, int particle, int q_num, int d) {

    Walker* diff_walker = new Walker(n_p, dim, false);
    diff_walker->r = walker->r;
    diff_walker->r(particle, d) += h;
    diff_walker->calc_r_i2(particle);
    set_qnum_indie_terms(diff_walker, particle);
    double phi_pluss = basis_functions[q_num]->eval(diff_walker, particle);

    diff_walker->r(particle, d) -= 2 * h;
    diff_walker->calc_r_i2(particle);
    set_qnum_indie_terms(diff_walker, particle);
    double phi_minus = basis_functions[q_num]->eval(diff_walker, particle);

    delete diff_walker;

    return (phi_pluss - phi_minus) / (2 * h);

}

double Orbitals::num_ddiff(const Walker* walker, int particle, int q_num) {
    Walker* diff_walker = new Walker(n_p, dim, false);
    diff_walker->r = walker->r;

    double ddiff = 0;
    double phi_orig = walker->phi(particle, q_num);
    for (int k = 0; k < dim; k++) {

        diff_walker->r(particle, k) += h;
        diff_walker->calc_r_i2(particle);
        set_qnum_indie_terms(diff_walker, particle);
        ddiff += basis_functions[q_num]->eval(diff_walker, particle);

        diff_walker->r(particle, k) -= 2 * h;
        diff_walker->calc_r_i2(particle);
        set_qnum_indie_terms(diff_walker, particle);
        ddiff += basis_functions[q_num]->eval(diff_walker, particle);

        diff_walker->r(particle, k) += h;
        diff_walker->r2(particle) = walker->r2(particle);
    }

    //contribution from phi(r, q_num) is constant in loop and already tabulated.
    ddiff = (ddiff - 2 * dim * phi_orig) * h2;

    delete diff_walker;

    return ddiff;
}

//double Orbitals::get_variational_derivative(const Walker* walker, int n) {
//
//    Walker* diff_walker = new Walker(n_p, dim);
//    diff_walker->r = walker->r;
//    diff_walker->r2 = walker->r2;
//
//    double a = get_parameter(n);
//    set_parameter(a + h, 0);
//
//    qmc->get_sampling_ptr()->set_trial_states(diff_walker);
//    double phip = qmc->get_system_ptr()->get_spatial_wf(diff_walker);
//
//    set_parameter(a - h, 0);
//    qmc->get_sampling_ptr()->set_trial_states(diff_walker);
//    double phim = qmc->get_system_ptr()->get_spatial_wf(diff_walker);
//
//    set_parameter(a, 0);
//
//    return (phip - phim) / (2 * h * qmc->get_system_ptr()->get_spatial_wf(walker));
//}

double Orbitals::get_variational_derivative(Walker* walker, double alpha) {

    double dell_alpha = 0;

    //Fermions only
    for (int i = 0; i < n_p; i++) {
        for (int j = 0; j < n2; j++) {
            dell_alpha += walker->phi(i, j) * get_dell_alpha_phi(walker, i, j, alpha) * walker->inv(j, i);
        }
    }

    return dell_alpha;

}

double Orbitals::get_dell_alpha_phi(Walker* walker, int p, int q_num, int n){
    Walker* diff_walker = new Walker(n_p, dim);

    double alpha_orig = get_parameter(n);

    double phi_mid = walker->phi(p, q_num);

    diff_walker->r.row(p) = walker->r.row(p);
    diff_walker->r2(p) = walker->r2(p);

    set_parameter(alpha_orig + h, n);
    set_qnum_indie_terms(diff_walker, p);
    double phi_plus = phi(diff_walker, p, q_num);

    set_parameter(alpha_orig - h, n);
    set_qnum_indie_terms(diff_walker, p);
    double phi_minus = phi(diff_walker, p, q_num);

    delete diff_walker;
    set_parameter(alpha_orig, n);
    set_qnum_indie_terms(walker, p);

    
    
    return (phi_plus - phi_minus)/(2*h*phi_mid);


}

double Orbitals::phi(const Walker* walker, int particle, int q_num) {
    return basis_functions[q_num]->eval(walker, particle);
}

double Orbitals::del_phi(const Walker* walker, int particle, int q_num, int d) {
//        testDell(walker, particle, q_num, d);
    return dell_basis_functions[d][q_num]->eval(walker, particle);
}

double Orbitals::lapl_phi(const Walker* walker, int particle, int q_num) {
//        testLaplace(walker, particle, q_num);
    return lapl_basis_functions[q_num]->eval(walker, particle);
}

void Orbitals::testDell(const Walker* walker, int particle, int q_num, int d) {
    double cf = dell_basis_functions[d][q_num]->eval(walker, particle);
    double num = num_diff(walker, particle, q_num, d);
    double diff = fabs(cf / num - 1);
    if (diff > 1e-2) {
        std::cout << "dell " << particle << "  " << q_num << "  " << d;
        std::cout << "   " << diff << std::endl;
    }
}

void Orbitals::testLaplace(const Walker* walker, int particle, int q_num) {
    double cf = lapl_basis_functions[q_num]->eval(walker, particle);
    double num = num_ddiff(walker, particle, q_num);
    double diff = fabs(cf / num - 1);
    if (diff > 1e-2) {
        std::cout << "lapl " << particle << "  " << q_num;
        std::cout << "      " << diff << std::endl;

    }
}
