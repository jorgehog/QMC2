/* 
 * File:   VMC.cpp
 * Author: jorgmeister
 * 
 * Created on October 12, 2012, 2:43 PM
 */

#include "../../QMCheaders.h"

VMC::VMC(GeneralParams & gP, VMCparams & vP, SystemObjects & sO, ParParams & pp)
: QMC(gP, vP.n_c, sO, pp) {

    name = "vmc";

    original_walker = new Walker(n_p, dim);

    pop_tresh = vP.pop_tresh;
    last_walker = 0;

    vmc_E = 0;
    thermalization = n_c / 10 * (n_c < 1e6) + 1e5 * (n_c >= 1e6);

    sampling->set_dt(vP.dt);

    set_trial_positions();
    copy_walker(original_walker, trial_walker);
    
}

VMC::~VMC() {
    delete trial_walker;
    delete [] original_walkers;
}

void VMC::set_trial_positions() {
    sampling->set_trial_pos(original_walker);

}

void VMC::store_walker() {
    if (cycle % pop_tresh == 0){
        copy_walker(original_walker, original_walkers[last_walker]);
        original_walkers[last_walker]->set_E(local_E);
        last_walker++;
    }
}

//void VMC::initialize() {
//
//    sampling->set_trial_pos(original_walker);
//
//}

void VMC::run_method() {

    //    initialize();

    for (cycle = 1; cycle <= thermalization; cycle++) {
        diffuse_walker(original_walker, trial_walker);

    }

    for (cycle = 1; cycle <= n_c; cycle++) {
        
        diffuse_walker(original_walker, trial_walker);

        calculate_energy_necessities(original_walker);
        local_E = calculate_local_energy(original_walker);

        vmc_E += local_E;

        dump_output();
        error_estimator->update_data(local_E);
        store_walker();
        
    }


    dump_distribution();
    node_comm();
    scale_values();
    output();
    finalize_output();
    estimate_error();
}

void VMC::output() {
    using namespace std;

    s << "VMC energy: " << get_energy() << endl;
    s << "Acceptance ratio: " << get_accepted_ratio(n_p * n_w * (thermalization + n_c)) << endl;

    std_out->cout(s);
}

void VMC::node_comm() {
#ifdef MPI_ON
    MPI_Allreduce(MPI_IN_PLACE, &vmc_E, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
}


