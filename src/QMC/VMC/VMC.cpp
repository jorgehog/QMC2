/* 
 * File:   VMC.cpp
 * Author: jorgmeister
 * 
 * Created on October 12, 2012, 2:43 PM
 */

#include "../../QMCheaders.h"

VMC::VMC(GeneralParams & gP, VMCparams & vP, SystemObjects & sO, ParParams & pp, int n_w)
: QMC(gP, vP.n_c, sO, pp, n_w) {

    name = "vmc";
    original_walker = new Walker(n_p, dim);

    dist_tresh = 50;
    pop_tresh = n_c / n_w;
    last_walker = 0;

    vmc_E = 0;
    thermalization = n_c / 10 * (n_c < 1e6) + 1e5 * (n_c >= 1e6);

    sampling->set_dt(vP.dt);

}

VMC::~VMC() {
    delete trial_walker;
    delete [] original_walkers;
}

void VMC::set_trial_positions() {
    sampling->set_trial_pos(original_walker);

}

void VMC::save_distribution() {
//    using namespace std;
    if (cycle % dist_tresh == 0) {
//        cout << "-----preDist------" << endl;
//                cout << dist << endl; 
//        cout << "-----thisR--------"<< endl;
//                cout << original_walker->r << endl;
        dist.insert_rows(dist.n_rows, original_walker->r);


//        arma::mat rlol = arma::zeros<arma::mat > (2, 2);
//        rlol(0, 0) = 1;
//        rlol(0, 1) = -1;
//        rlol(1, 0) = 1;
//        rlol(1, 1) = -1;
//        if (cycle == dist_tresh) {
//            std::cout << rlol << std::endl;
//        }
//        dist.insert_rows(dist.n_rows, rlol);

//        cout << "-----newDist----" << endl;
//        cout << dist << endl;
//        cout << "-----------------" << endl;
//        sleep(10);
    }
}

void VMC::store_walkers() {

    //store for DMC
    if (cycle % pop_tresh == 0) {
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

        dump_output();
        error_estimator->update_data(local_E);
        store_walkers();

    }

    node_comm();
    scale_values();
    output();
    finalize_output();
    estimate_error();
}

void VMC::output() {
    using namespace std;

    s << "VMC energy: " << get_energy() << endl;
    std_out->cout(s);
    get_accepted_ratio();

}

void VMC::node_comm() {
#ifdef MPI_ON
    MPI_Allreduce(MPI_IN_PLACE, &vmc_E, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
}


