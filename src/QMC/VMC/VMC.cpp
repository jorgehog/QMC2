/* 
 * File:   VMC.cpp
 * Author: jorgmeister
 * 
 * Created on October 12, 2012, 2:43 PM
 */

#include "../../QMCheaders.h"


VMC::VMC(GeneralParams & gP, VMCparams & vP, SystemObjects & sO)
: QMC(gP.n_p, gP.dim, vP.n_c, sO.sample_method, sO.SYSTEM, sO.kinetics, sO.jastrow) {

    vmc_E = 0;
    E2 = 0;

    original_walker = new Walker(n_p, dim);
    trial_walker = new Walker(n_p, dim);

    this->thermalization = n_c / 10 * (n_c < 1e6) + 1e5 * (n_c >= 1e6);

}

void VMC::initialize() {
    jastrow->initialize();

    sampling->set_trial_pos(original_walker);

    copy_walker(original_walker, trial_walker);
}

bool VMC::move_autherized(double A) {
    return metropolis_test(A);
}

void VMC::calculate_energy(Walker* walker) {

    local_E = calculate_local_energy(walker);

    vmc_E += local_E;
    E2 += local_E*local_E;
}

void VMC::scale_values() {
    vmc_E /= n_c;
    E2 /= n_c;
}

void VMC::run_method() {

    initialize();
    
    for (cycle = 1; cycle <= thermalization; cycle++) {
        diffuse_walker(original_walker, trial_walker);
    }

    for (cycle = 1; cycle <= n_c; cycle++) {

        diffuse_walker(original_walker, trial_walker);

        calculate_energy_necessities(original_walker);
        calculate_energy(original_walker);

        dump_output();

    }

    scale_values();
    user_output();
    finalize_output();

}

void VMC::user_output() const {
    std::cout << "VMC energy: " << get_energy() << std::endl;
    std::cout << "VMC variance: " << get_var() << std::endl;
    std::cout << "Acceptance ratio: " << get_accepted_ratio(n_p * (thermalization + n_c)) << endl;
}

double VMC::get_var() const {
    return E2 - vmc_E*vmc_E;
}

double VMC::get_e2() const {
    return E2;
}

void VMC::set_e(double E) {
    vmc_E = E;
}

void VMC::set_e2(double E2) {
    this->E2 = E2;
}

double VMC::get_energy() const {
    return vmc_E;
}
