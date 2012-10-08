/* 
 * File:   QMC.cpp
 * Author: jorgehog
 * 
 * Created on 30. mars 2012, 17:42
 */

#include "../QMCheaders.h"

QMC::QMC(int n_p, int dim, int n_c,
        Sampling *sampling,
        System *system,
        Kinetics *kinetics,
        Jastrow *jastrow) {

    this->n_p = n_p;
    this->dim = dim;
    this->n_c = n_c;
    this->n2 = n_p / 2;

    this->jastrow = jastrow;
    this->sampling = sampling;
    this->system = system;
    this->kinetics = kinetics;

    this->sampling->set_qmc_ptr(this);
    this->kinetics->set_qmc_ptr(this);

    this->accepted = 0;

}

QMC::QMC() {

}

void QMC::add_output(OutputHandler* output_handler) {
    output_handler->set_qmc_ptr(this);
    this->output_handler.push_back(output_handler);
}

void QMC::dump_output() {

    for (std::vector<OutputHandler*>::iterator output_obj = output_handler.begin(); output_obj != output_handler.end(); ++output_obj) {
        (*output_obj)->dump();
    }

}

void QMC::finalize_output() {

    for (std::vector<OutputHandler*>::iterator output_obj = output_handler.begin(); output_obj != output_handler.end(); ++output_obj) {
        (*output_obj)->finalize();
    }

}

void QMC::get_gradients(Walker* walker, int particle) const {
    jastrow->get_grad(walker);
    system->get_spatial_grad(walker, particle);
}

void QMC::get_gradients(Walker* walker) const {
    jastrow->get_grad(walker);
    system->get_spatial_grad(walker, 0);
    system->get_spatial_grad(walker, n2);
}

void QMC::get_laplsum(Walker* walker) const {
    walker->lapl_sum = system->get_spatial_lapl_sum(walker) + jastrow->get_lapl_sum(walker);
}

void QMC::get_wf_value(Walker* walker) const {
    walker->value = system->get_spatial_wf(walker) * jastrow->get_val(walker);
}

void QMC::update_pos(const Walker* walker_pre, Walker* walker_post, int particle) const {

    for (int j = 0; j < dim; j++) {
        walker_post->r(particle, j) = walker_pre->r(particle, j)
                + sampling->get_new_pos(walker_pre, particle, j);
    }

    for (int j = 0; j < n_p; j++) {
        if (j != particle) {
            walker_post->r_rel(particle, j) = walker_post->r_rel(j, particle)
                    = walker_post->abs_relative(particle, j);
        }
    }

    walker_post->calc_r_i2(particle);

    sampling->update_necessities(walker_pre, walker_post, particle);

}

double QMC::get_acceptance_ratio(const Walker* walker_pre, const Walker* walker_post, int particle) const {
    double spatial_jast = sampling->get_spatial_ratio(walker_post, walker_pre, particle);
    double G = sampling->get_g_ratio(walker_post, walker_pre);

    return spatial_jast * spatial_jast * G;
}

void QMC::calculate_energy_necessities(Walker* walker) const {
    kinetics->calculate_energy_necessities(walker);
}

bool QMC::metropolis_test(double A) {
    double r = sampling->call_RNG();

    if (r <= A) {
        accepted++;
        return true;

    } else {
        return false;
    }
}

void QMC::update_walker(Walker* walker_pre, const Walker* walker_post, int particle) const {

    for (int i = 0; i < dim; i++) {
        walker_pre->r(particle, i) = walker_post->r(particle, i);
    }

    for (int i = 0; i < n_p; i++) {
        walker_pre->r_rel(i, particle) = walker_pre->r_rel(particle, i) = walker_post->r_rel(i, particle);
    }

    walker_pre->r2[particle] = walker_post->r2[particle];

    sampling->update_walker(walker_pre, walker_post, particle);
}

void QMC::reset_walker(const Walker* walker_pre, Walker* walker_post, int particle) const {
    for (int i = 0; i < dim; i++) {
        walker_post->r(particle, i) = walker_pre->r(particle, i);
    }

    for (int i = 0; i < n_p; i++) {
        walker_post->r_rel(i, particle) = walker_post->r_rel(particle, i) = walker_pre->r_rel(i, particle);
    }

    walker_post->r2[particle] = walker_pre->r2[particle];

    sampling->reset_walker(walker_pre, walker_post, particle);
}

void QMC::diffuse_walker(Walker* original, Walker* trial) {
    for (int particle = 0; particle < n_p; particle++) {

        update_pos(original, trial, particle);

        double A = get_acceptance_ratio(original, trial, particle);

        if (move_autherized(A)) {
            update_walker(original, trial, particle);
        } else {
            reset_walker(original, trial, particle);
        }

    }
}

void QMC::copy_walker(const Walker* parent, Walker* child) const {

    child->r2 = parent->r2;
    child->r = parent->r;
    child->r_rel = parent->r_rel;

    child->ressurect();
    child->set_E(parent->get_E());

    sampling->copy_walker(parent, child);

}

double QMC::calculate_local_energy(Walker* walker) const {
    return kinetics->get_KE(walker) + system->get_potential_energy(walker);
}

/*
 
 
 VMC
 
 
 */



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

/*
 
 
DMC 
 
 
 */



DMC::DMC(GeneralParams & gP, DMCparams & dP, SystemObjects & sO)
: QMC(gP.n_p, gP.dim, dP.n_c, sO.sample_method, sO.SYSTEM, sO.kinetics, sO.jastrow) {

    this->dist_from_file = dP.dist_in;

    this->block_size = dP.n_b;
    this->n_w = dP.n_w;
    this->thermalization = dP.therm;
    this->E_T = dP.E_T;

    K = 5;
    int max_walkers = K * n_w;

    original_walkers = new Walker*[max_walkers];
    trial_walker = new Walker(n_p, dim);

}

bool DMC::move_autherized(double A) {
    return metropolis_test(A)&(A > 0);
}

void DMC::initialize() {

    jastrow->initialize();

    //Initializing active walkers
    for (int k = 0; k < n_w; k++) {
        original_walkers[k] = new Walker(n_p, dim);
    }

    //Seting trial position of active walkers
    if (dist_from_file) {

        ifstream dist;
        string path = "/home/jorgmeister/scratch/";
        string name = "dist_out.dat";
        dist.open((path + name).c_str());

        for (int k = 0; k < n_w; k++) {
            sampling->set_trial_pos(original_walkers[k], true, &dist);
        }

        dist.close();

    } else {

        for (int k = 0; k < n_w; k++) {
            while (original_walkers[k]->is_singular()) {
                sampling->set_trial_pos(original_walkers[k]);
            }
        }

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

void DMC::user_output() const {
    printf("dmcE: %1.5f| Nw: %4d| %1.5f%%", dmc_E / cycle, n_w,
            (double) cycle / n_c * 100);
    cout << endl;
}

void DMC::Evolve_walker(int k, double GB) {

    int branch_mean = int(GB + sampling->call_RNG());
    double dE = (original_walkers[k]->get_E() - E_T);
    dE = dE*dE;

    if (branch_mean == 0 || dE > 1. / sampling->get_dt()) {
        original_walkers[k]->kill();
        //        cout << "died" << k << endl;

    } else {

        for (int n = 1; n < branch_mean; n++) {
            copy_walker(original_walkers[k], original_walkers[n_w]);
            //            cout << "spawned" << k << " to " << n_w << endl;
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
        //        if ((abs(local_E_prev - 20.16) > 0.1) && n_b != 1) {
        //            cout << local_E_prev << endl;
        //        }
        diffuse_walker(original_walkers[k], trial_walker);
        //        cout << k << " "<<n_w << " "<<cycle <<" "<<b<<endl;
        //        if ((k == 981) && (n_w == 982) && (cycle == 7) && (b == 0)) {
        //            original_walkers[k]->print("ETTER");
        //        }
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

        user_output();

    }


    dmc_E = 0;
    for (cycle = 1; cycle <= n_c; cycle++) {

        reset_parameters();

        for (int k = 0; k < n_w_last; k++) {
            iterate_walker(k, block_size);
        }

        bury_the_dead();
        update_energies();

        user_output();
        dump_output();

    }

    finalize_output();
}

void DMC::bury_the_dead() {

    int newborn = n_w - n_w_last;
    int last_alive = n_w - 1;
    int i = 0;
    int k = 0;
    bool cond = false;
    if (cond) cout << "--init-- NEWBORN: " << newborn << " DEAD: " << deaths << " last alive: " << last_alive << endl;
    while (k < newborn && i < n_w_last) {
        if (original_walkers[i]->is_dead()) {
            //            cout << "--spawned replace-- walker[" << i << "] is dead. replacing with" << last_alive << endl;
            copy_walker(original_walkers[last_alive], original_walkers[i]);
            delete original_walkers[last_alive];
            original_walkers[last_alive] = new Walker(n_p, dim, false);
            last_alive--;
            k++;
        }
        i++;
    }
    if (cond) cout << "--spawned replace--" << k << " replacements made. Last alive is now " << last_alive << endl;


    if (deaths > newborn) {
        int difference = deaths - newborn;
        if (cond) cout << "--remaining dead--we must now shuffle " << difference << " walkers" << endl;
        int i = 0;
        int first_dead = 0;

        //we have to delete [difference] spots to compress array.
        while (i != difference) {

            //Finds last living walker and deletes any dead walkers at array end
            while (original_walkers[last_alive]->is_dead()) {
                //                cout << "--remaining dead--dead walker at endpoint: " << last_alive << endl;
                delete original_walkers[last_alive];
                original_walkers[last_alive] = new Walker(n_p, dim, false);
                i++;
                last_alive--;
            }
            if (cond) cout << "--remaining dead--last alive walker found: " << last_alive << endl;
            //Find the first dead walker
            while (original_walkers[first_dead]->is_alive()) {
                first_dead++;
            }
            if (cond) cout << "--remaining dead--first dead found: " << first_dead << endl;
            //If there is a dead walker earlier in the array, the last living
            //takes the spot.
            if (first_dead < last_alive) {
                if (cond) cout << "--remaining dead--walker[" << first_dead << "] is dead. replacing with" << last_alive << endl;
                copy_walker(original_walkers[last_alive], original_walkers[first_dead]);
                delete original_walkers[last_alive];
                original_walkers[last_alive] = new Walker(n_p, dim, false);

                i++;
                last_alive--;
                first_dead++;
            }
        }
        if (cond) cout << "--remaining dead--" << i << " replacements made. " << endl;
    }


    n_w = n_w - deaths;
    if (cond) cout << "--fin--" << "Last alive is now " << last_alive << " NW now: " << n_w << endl;

}