/* 
 * File:   VMC.cpp
 * Author: jorgmeister
 * 
 * Created on October 12, 2012, 2:43 PM
 */

#include "../../QMCheaders.h"

VMC::VMC(GeneralParams & gP, VMCparams & vP, SystemObjects & sO)
: QMC(gP.n_p, gP.dim, vP.n_c, sO.sample_method, sO.SYSTEM, sO.jastrow) {

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
        error_estimator->update_data(local_E);

    }

    scale_values();
    user_output();
    finalize_output();
    estimate_error();
}

void VMC::user_output() const {
    using namespace std;

    cout << "VMC energy: " << get_energy() << endl;
    cout << "VMC variance: " << get_var() << endl;
    cout << "Acceptance ratio: " << get_accepted_ratio(n_p * (thermalization + n_c)) << endl;
}


