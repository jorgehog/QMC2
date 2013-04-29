/* 
 * File:   VMC.cpp
 * Author: jorgmeister
 * 
 * Created on October 12, 2012, 2:43 PM
 */

#include "../../QMCheaders.h"

VMC::VMC(GeneralParams & gP, VMCparams & vP, SystemObjects & sO, ParParams & pp, int n_w, bool dist_out)
: QMC(gP, vP.n_c, sO, pp, n_w) {

    output_tresh = n_c / 100;

    if (output_tresh == 0) {
        output_tresh = 1;
    }

    original_walker = new Walker(n_p, dim);

    dist_tresh = 25;
    pop_tresh = n_c / n_w;
    offset = n_c - n_w*pop_tresh;

    last_walker = 0;

    vmc_E = 0;
    thermalization = n_c / 10 * (n_c < 1e6) + 1e5 * (n_c >= 1e6);

    sampling->set_dt(vP.dt);

    if (dist_out) {
        dist = arma::zeros(n_c * n_p / dist_tresh, dim);
    }

}

VMC::~VMC() {
    delete trial_walker;
    delete [] original_walkers;
}

void VMC::set_trial_positions() {
    sampling->set_trial_pos(original_walker);

}

void VMC::save_distribution() {
    using namespace arma;

    if (cycle % dist_tresh == 0) {
        dist(span(last_inserted, last_inserted + n_p - 1), span()) = original_walker->r;
        last_inserted += n_p;
    }
}

void VMC::store_walkers() {

    if ((cycle > offset) && (cycle % pop_tresh == 0)) {
        copy_walker(original_walker, original_walkers[last_walker]);
        original_walkers[last_walker]->set_E(local_E);
        last_walker++;
    }

}

void VMC::run_method() {

    set_trial_positions();
    copy_walker(original_walker, trial_walker);

    for (cycle = 1; cycle <= thermalization; cycle++) {
        diffuse_walker(original_walker, trial_walker);
    }

    for (cycle = 1; cycle <= n_c; cycle++) {

        diffuse_walker(original_walker, trial_walker);

        calculate_energy_necessities(original_walker);
        local_E = calculate_local_energy(original_walker);

        vmc_E += local_E;


        output();
        update_subsamples();
        dump_output();
        error_estimator->update_data(local_E);
        store_walkers();
    }

    node_comm();
    scale_values();

    output();
    dump_subsamples();
    get_accepted_ratio();

    finalize_output();
    estimate_error();
}

void VMC::output() {
    using namespace std;

    if (is_master) {
        if (cycle > n_c) {
            cout << endl;
            cout << setprecision(6) << fixed;
            cout << "Final VMC energy: " << vmc_E << endl;
        } else if (cycle % output_tresh == 0) {
            cout << "\rVMC energy: " << vmc_E / cycle << "  " << (double) cycle / n_c * 100 << "%";
            cout.flush();
        }

    }

}

void VMC::node_comm() {
#ifdef MPI_ON
    MPI_Allreduce(MPI_IN_PLACE, &vmc_E, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
}


