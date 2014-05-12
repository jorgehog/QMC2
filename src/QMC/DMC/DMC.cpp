#include "DMC.h"

#include <iomanip>

#include "../../structs.h"

#include "../../System/System.h"
#include "../../Sampling/Sampling.h"
#include "../../Walker/Walker.h"
#include "../../OutputHandler/stdoutDMC/stdoutDMC.h"
#include "../../OutputHandler/Distribution/Distribution.h"
#include "../../ErrorEstimator/ErrorEstimator.h"
#include "../VMC/VMC.h"
#include "../../Orbitals/Orbitals.h"


using namespace QMC2;


DMC::DMC(GeneralParams & gP, DMCparams & dP, SystemObjects & sO, ParParams & pp, bool silent)
: QMC(gP, dP.n_c, sO, pp, silent, dP.dt, dP.n_w, K) {

    if (is_master) DMCout = new stdoutDMC(this, gP.runpath);

    std::stringstream name;
    name << sO.SP_basis->getName() << gP.n_p << "c" << gP.systemConstant << "dmc";
    distribution = new Distribution(pp, gP.runpath, name.str(), silent);

    dist_tresh = 25;

    dist_size = n_w * n_p * n_c / dist_tresh;
    dist = arma::zeros(dist_size, dim);

    block_size = dP.n_b;
    thermalization = dP.therm;

    thermalized = false;
    force_comm = false;

    dmc_E = 0;
    dmc_E_unscaled = 0;

    initializedFromVMC = false;

    if (parallel) {
        n_w_list = arma::zeros<arma::uvec > (n_nodes);
    }


}

void DMC::initFromVMC(VMC *vmc)
{
    E_T = vmc->get_energy();

    for (unsigned int i = 0; i < n_w; i++) {
        copy_walker(vmc->original_walkers[i], original_walkers[i]);
        //delete vmc->original_walkers[i];
    }
    //delete [] vmc->original_walkers;

    initializedFromVMC = true;

}

void DMC::set_trial_positions() {

    if (initializedFromVMC) return;

    double tmpDt = sampling->get_dt();
    sampling->set_dt(0.5);
    for (unsigned int k = 0; k < n_w; k++) {
        sampling->set_trial_pos(original_walkers[k]);
    }
    sampling->set_dt(tmpDt);

    //Calculating and storing energies of active walkers
    E_T = 0;
    for (unsigned int k = 0; k < n_w; k++) {
        calculate_energy_necessities(original_walkers[k]);
        double El = calculate_local_energy(original_walkers[k]);

        original_walkers[k]->set_E(El);
        E_T += El;
    }

    E_T /= n_w;

}

void DMC::output() {

    using namespace std;

    stringstream s;

    if (is_master) {

        s << setprecision(6) << fixed;
        s << "dmcE: " << dmc_E << " | E_T: " << E_T;

        s << " | Nw: ";
        s << setfill(' ') << setw(5);
        s << n_w_tot << " | ";

        s << setprecision(1) << fixed << setfill(' ') << setw(5);
        s << (double) cycle / n_c * 100 << "%";

        if ((cycle == n_c) && thermalized) {
            cout << "\r" << s.str() << endl;
        } else {
            cout << "\r" << s.str();
            cout.flush();

#ifdef LINED_OUTPUT
                cout << endl;
#endif

        }

    }

}

void DMC::Evolve_walker(int k, double GB) {

    int branch_mean = int(GB + sampling->call_RNG());

    if (branch_mean == 0) {
        original_walkers[k]->kill();
    } else {

        for (int n = 1; n < branch_mean; n++) {
            copy_walker(original_walkers[k], original_walkers[n_w]);
            n_w++;
        }

        E += GB * local_E;
        samples++;

        if (thermalized) {
            update_subsamples(GB);
        }

    }
}

void DMC::update_energies() {

    node_comm();

    E_T = E / samples;
    dmc_E_unscaled += E_T;
    dmc_E = dmc_E_unscaled / cycle;

    push_subsamples();


}

void DMC::iterate_walker(int k) {

    copy_walker(original_walkers[k], trial_walker);

    for (unsigned int b = 0; b < block_size; b++) {

        double local_E_prev = original_walkers[k]->get_E();

        diffuse_walker(original_walkers[k], trial_walker);

        if (thermalized) update_samplers(original_walkers[k]);

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

void DMC::run_method(bool initialize) {

    sampling->set_dt(dtOrig);

    if (initialize) {
        set_trial_positions();
    } else {
        reset_all();
    }


    for (cycle = 1; cycle <= thermalization; cycle++) {

        reset_parameters();

        for (unsigned int k = 0; k < n_w_last; k++) {
            iterate_walker(k);
        }

        bury_the_dead();

        update_energies();

        normalize_population();

        output();


    }

    normalize_population();

    thermalized = true;
    dmc_E = dmc_E_unscaled = 0;

    for (cycle = 1; cycle <= n_c; cycle++) {

        reset_parameters();

        for (unsigned int k = 0; k < n_w_last; k++) {
            iterate_walker(k);
        }

        error_estimator->update_data(E / samples);

        bury_the_dead();
        update_energies();

        normalize_population();

        output();
        save_distribution();
        if (is_master) DMCout->dump();

    }

    if (is_master) DMCout->finalize();

//    free_walkers();
//    dump_subsamples(true); //This should be called from the outside if wanted.
    estimate_error();
    finalize_distribution();
    get_accepted_ratio();
    clean();
}

void DMC::save_distribution() {
    using namespace arma;

    if (cycle % dist_tresh == 0) {

        if ((last_inserted + n_p * n_w) >= dist.n_rows) {
            dist.resize(dist.n_rows + n_w * n_p * 100, dim);
        }

        for (unsigned int i = 0; i < n_w; i++) {
            dist(span(last_inserted, last_inserted + n_p - 1), span()) = original_walkers[i]->r;
            last_inserted += n_p;
        }

    }
}

void DMC::bury_the_dead() {

    unsigned int newborn = n_w - n_w_last;
    unsigned int last_alive = n_w - 1;
    unsigned int i = 0;
    unsigned int k = 0;

    while (k < newborn && i < n_w_last) {
        if (original_walkers[i]->is_dead()) {

            copy_walker(original_walkers[last_alive], original_walkers[i]);
            original_walkers[last_alive]->kill();
            //            delete original_walkers[last_alive];
            //            original_walkers[last_alive] = new Walker(n_p, dim, false);
            last_alive--;
            k++;
        }
        i++;
    }


    if (deaths > newborn) {

        unsigned int difference = deaths - newborn;
        unsigned int i = 0;
        unsigned int first_dead = 0;

        //we have to delete [difference] spots to compress array.
        while (i != difference) {

            //Finds last living walker and deletes any dead walkers at array end
            while (original_walkers[last_alive]->is_dead()) {

                //                delete original_walkers[last_alive];
                //                original_walkers[last_alive] = new Walker(n_p, dim, false);
                original_walkers[last_alive]->kill();
                i++;

                //if all the walkers are dead
                if (last_alive == 0) {
                    n_w = 0;
                    return;
                }

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
                //                delete original_walkers[last_alive];
                //                original_walkers[last_alive] = new Walker(n_p, dim, false);
                original_walkers[last_alive]->kill();

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

    if (parallel) {

        MPI_Allreduce(MPI_IN_PLACE, &E, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &samples, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allgather(&n_w, 1, MPI_INT, n_w_list.memptr(), 1, MPI_INT, MPI_COMM_WORLD);

        for (unsigned int i = 0; i < n_nodes; i++) {
            if (n_w_list(i) == 0) {
                force_comm = true;
            }
        }

        n_w_tot = arma::accu(n_w_list);

    } else {
        n_w_tot = n_w;
    }
#else
    n_w_tot = n_w;
#endif

    if (n_w_tot == 0) {
        if (is_master) {
            std::cout << "All walkers dead. Exiting..." << std::endl;
            exit(1);
        }
    }

}

void DMC::switch_souls(unsigned int root, unsigned int root_id, unsigned  int dest, unsigned  int dest_id) {
    if (node == root) {
        original_walkers[root_id]->send_soul(dest);
        n_w--;
    } else if (node == dest) {
        original_walkers[dest_id]->recv_soul(root);
        n_w++;
    }
}

void DMC::normalize_population() {
#ifdef MPI_ON
    using namespace arma;

    if (!(cycle % check_thresh == 0) && !force_comm) return;

    force_comm = false;

    unsigned int avg = n_w_tot / n_nodes;

    umat swap_map = zeros<umat > (n_nodes, n_nodes); //root x (recieve_count @ index dest)
    uvec snw = sort_index(n_w_list, 1); //enables us to index n_w as decreasing

    //Start iterating sending from highest to lowest. When root or dest reaches the average
    //value they are shifted to the second highest/lowest. Process continues untill root
    //reaches dest.
    unsigned int root = 0;
    unsigned int dest = n_nodes - 1;
    while (root < dest) {
        if (n_w_list(snw(root)) > avg) {
            if (n_w_list(snw(dest)) < avg) {
                swap_map(snw(root), snw(dest))++;
                n_w_list(snw(root))--;
                n_w_list(snw(dest))++;
            } else {
                dest--;
            }
        } else {
            root++;
        }
    }

    //At this point we have at most a unbalance of n_nodes-1 at one node.
    //optimal in this case is a unbalance of 1 on n_nodes-1 nodes.
    //To do this shuffeling we send walkers from the highest to the lowest untill
    //the highest reaches avg+1 walkers. The lowest is shifted every time it recieves a 
    //walker.

    snw = sort_index(n_w_list, 1);
    root = 0;
    dest = n_nodes - 1;
    while (root < dest) {
        if (n_w_list(snw(root)) > avg + 1) {

            swap_map(snw(root), snw(dest))++;
            n_w_list(snw(root))--;
            n_w_list(snw(dest))++;

            dest--;

        } else {
            root++;
        }
    }

    uvec test = sum(swap_map, 1);
    if (test.max() < sendcount_thresh) {
        test.reset();
        swap_map.reset();
        return;
    }


    for (unsigned int root = 0; root < n_nodes; root++) {
        for (unsigned int dest = 0; dest < n_nodes; dest++) {
            if (swap_map(root, dest) != 0) {

                for (unsigned int sendcount = 0; sendcount < swap_map(root, dest); sendcount++) {
                    switch_souls(root, n_w - 1, dest, n_w);
                }


            }
        }
    }

    test.reset();
    swap_map.reset();
    MPI_Barrier(MPI_COMM_WORLD);

#endif
}

void DMC::free_walkers() {

    delete trial_walker;
    for (unsigned int i = 0; i < n_w_size; i++) {
        delete original_walkers[i];
    }

    delete [] original_walkers;

}

void DMC::reset_all()
{

    thermalization = 0;

    if (is_master) DMCout->reset();

    QMC::reset_all();

}


bool DMC::move_autherized(double A) {
    return metropolis_test(A) & system->allow_transition();
}
