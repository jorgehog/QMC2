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
    using namespace std;
    initialize();

    //    cout << original_walker->r << endl;
    //    cout << original_walker->r2 << endl;
    //    
    //    for (int i = 0; i < n2; i++) {
    //        get_orbitals_ptr()->set_qnum_indie_terms(original_walker, 0);
    //        cout << "phi" << i << " " << get_orbitals_ptr()->phi(original_walker, 0, i) << endl;
    //        for (int k = 0; k < dim; k++) {
    //            cout << "dphi" << i << " "<<k <<" " << get_orbitals_ptr()->del_phi(original_walker, 0, i, k) << endl;
    //        }
    //        cout << "lphi" << i << " " << get_orbitals_ptr()->lapl_phi(original_walker, 0, i) << endl;
    //    }
    //    exit(1);

    //    cout << original_walker->inv << endl;
    //    cout << arma::accu(original_walker->inv) << endl;
    //    cout << original_walker->jast_grad << endl;
    //    cout << original_walker->spatial_grad << endl;
    //    exit(1);

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
    using namespace std;

    cout << "VMC energy: " << get_energy() << endl;
    cout << "VMC variance: " << get_var() << endl;
    cout << "Acceptance ratio: " << get_accepted_ratio(n_p * (thermalization + n_c)) << endl;
}


