/* 
 * File:   DMC.cpp
 * Author: jorgmeister
 * 
 * Created on October 12, 2012, 2:42 PM
 */

#include "../../QMCheaders.h"

DMC::DMC(GeneralParams & gP, DMCparams & dP, SystemObjects & sO, ParParams & pp)
: QMC(gP.n_p, gP.dim, dP.n_c, sO, pp) {

    this->dist_from_file = dP.dist_in;
    this->dist_in_path = dP.dist_in_path;

    this->block_size = dP.n_b;
    this->n_w = dP.n_w;
    this->thermalization = dP.therm;
    this->E_T = dP.E_T;

    K = 5;
    int max_walkers = K * n_w;

    original_walkers = new Walker*[max_walkers];
    trial_walker = new Walker(n_p, dim);

}

void DMC::initialize() {
    jastrow->initialize();

    //Initializing active walkers
    for (int k = 0; k < n_w; k++) {
        original_walkers[k] = new Walker(n_p, dim);
    }

    //Seting trial position of active walkers
    if (dist_from_file) {


        for (int k = 0; k < n_w; k++) {
            s << dist_in_path << "walker_positions/dist_out" << node << "_" << k << ".arma";

            original_walkers[k]->r.load(s.str());
            sampling->set_trial_pos(original_walkers[k], false);
            

            s.str(std::string());
        }

    } else {
        double tmpDt = sampling->get_dt();
        sampling->set_dt(0.5);
        for (int k = 0; k < n_w; k++) {
            sampling->set_trial_pos(original_walkers[k]);
        }
        sampling->set_dt(tmpDt);

    }

    //Calculating and storing energies of active walkers
    for (int k = 0; k < n_w; k++) {
        calculate_energy_necessities(original_walkers[k]);
        double El = calculate_local_energy(original_walkers[k]);

        original_walkers[k]->set_E(El);
    }

    if (E_T == 0) {
        for (int k = 0; k < n_w; k++) {
            E_T += original_walkers[k]->get_E();
        }
        E_T /= n_w;
    }

    //Creating unactive walker objects (note: 3. arg=false implies dead) 
    for (int k = n_w; k < K * n_w; k++) {
        original_walkers[k] = new Walker(n_p, dim, false);
    }

}

void DMC::output() {
   s << "dmcE:" << dmc_E / cycle << "| Nw: " << n_w << "| " << (double) cycle / n_c * 100 << "%";
   std_out->cout(s);
}

void DMC::Evolve_walker(int k, double GB) {

    int branch_mean = int(GB + sampling->call_RNG());
    double dE = (original_walkers[k]->get_E() - E_T);
    dE = dE*dE;

    if (branch_mean == 0 || dE > 1. / sampling->get_dt()) {
        original_walkers[k]->kill();
    } else {

        for (int n = 1; n < branch_mean; n++) {
            copy_walker(original_walkers[k], original_walkers[n_w]);
            n_w++;
        }

        E += GB * local_E;
        samples++;

    }
}

void DMC::update_energies() {
    dmc_E += E / samples;
    E_T = dmc_E / cycle;
}

void DMC::iterate_walker(int k, int n_b) {

    copy_walker(original_walkers[k], trial_walker);

    for (int b = 0; b < n_b; b++) {

        double local_E_prev = original_walkers[k]->get_E();

        diffuse_walker(original_walkers[k], trial_walker);

        calculate_energy_necessities(original_walkers[k]);
        local_E = calculate_local_energy(original_walkers[k]);
        original_walkers[k]->set_E(local_E);

        double GB = sampling->get_branching_Gfunc(local_E, local_E_prev, E_T);

        Evolve_walker(k, GB);

        if (original_walkers[k]->is_dead()) {
            deaths++;
            break;
        }

    }
}

void DMC::run_method() {

    initialize();

    dmc_E = 0;
    for (cycle = 1; cycle <= thermalization; cycle++) {

        reset_parameters();

        for (int k = 0; k < n_w_last; k++) {
            iterate_walker(k, 1);
        }

        bury_the_dead();
        update_energies();

        output();

        node_comm();
    }
    
    cycle = 100;
    node_comm();

    dmc_E = 0;
    for (cycle = 1; cycle <= n_c; cycle++) {

        reset_parameters();

        for (int k = 0; k < n_w_last; k++) {
            iterate_walker(k, block_size);
        }

        bury_the_dead();
        update_energies();

        output();
        dump_output();

        error_estimator->update_data(E / samples);
        
        node_comm();

    }

    finalize_output();
    estimate_error();
}

void DMC::bury_the_dead() {

    int newborn = n_w - n_w_last;
    int last_alive = n_w - 1;
    int i = 0;
    int k = 0;


    while (k < newborn && i < n_w_last) {
        if (original_walkers[i]->is_dead()) {

            copy_walker(original_walkers[last_alive], original_walkers[i]);
            delete original_walkers[last_alive];
            original_walkers[last_alive] = new Walker(n_p, dim, false);
            last_alive--;
            k++;
        }
        i++;
    }


    if (deaths > newborn) {

        int difference = deaths - newborn;
        int i = 0;
        int first_dead = 0;

        //we have to delete [difference] spots to compress array.
        while (i != difference) {

            //Finds last living walker and deletes any dead walkers at array end
            while (original_walkers[last_alive]->is_dead()) {

                delete original_walkers[last_alive];
                original_walkers[last_alive] = new Walker(n_p, dim, false);
                i++;
                last_alive--;
            }

            //Find the first dead walker
            while (original_walkers[first_dead]->is_alive()) {
                first_dead++;
            }

            //If there is a dead walker earlier in the array, the last living
            //takes the spot.
            if (first_dead < last_alive) {

                copy_walker(original_walkers[last_alive], original_walkers[first_dead]);
                delete original_walkers[last_alive];
                original_walkers[last_alive] = new Walker(n_p, dim, false);

                i++;
                last_alive--;
                first_dead++;
            }
        }
    }

    n_w = n_w - deaths;

}

void DMC::node_comm() {
#ifdef MPI_ON
    if (cycle % 10 == 0){
        MPI_Allreduce(MPI_IN_PLACE, &E_T, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        if (is_master) {
            MPI_Reduce(MPI_IN_PLACE, &dmc_E, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            dmc_E/= n_nodes;
        } else {
            MPI_Reduce(&dmc_E, new double(), 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        }
        E_T /= n_nodes;
    } 
#endif
}