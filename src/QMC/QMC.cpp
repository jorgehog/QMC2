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
    this->thermalization = n_c / 10 * (n_c <= 1e6) + 1e5 * (n_c >= 1e6);

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

    update_necessities(walker_pre, walker_post, particle);

}

void QMC::update_necessities(const Walker* walker_pre, Walker* walker_post, int particle) const {
    sampling->update_necessities(walker_pre, walker_post, particle);
}

double QMC::get_acceptance_ratio(const Walker* walker_pre, const Walker* walker_post, int particle) const {
    double spatial_jast = sampling->get_spatial_ratio(walker_post, walker_pre, particle);
    double G = sampling->get_g_ratio(walker_post, walker_pre, particle);

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

        sampling->reset_control_parameters();
        update_pos(original, trial, particle);

        double A = get_acceptance_ratio(original, trial, particle);

        if (metropolis_test(A)) {
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

    sampling->copy_walker(parent, child);

}

double QMC::calculate_local_energy(Walker* walker) const {
    return kinetics->get_KE(walker) + system->get_potential_energy(walker);
}

/*
 
 
 VMC
 
 
 */
class VMC;

VMC::VMC(int n_p, int dim, int n_c,
        Sampling *sampling,
        System *system,
        Kinetics *kinetics,
        Jastrow *jastrow)
: QMC(n_p, dim, n_c, sampling, system, kinetics, jastrow) {

    vmc_E = 0;
    E2 = 0;

    original_walker = new Walker(n_p, dim);
    trial_walker = new Walker(n_p, dim);

}

void VMC::initialize() {
    jastrow->initialize();

    sampling->reset_control_parameters();
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

    }

    scale_values();

    finalize_output();

}

void VMC::output() const {
    std::cout << "VMC energy: " << get_energy() << std::endl;
    std::cout << "VMC variance: " << get_var() << std::endl;
    std::cout << "Acceptance ratio: " << get_accepted_ratio() << endl;
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


DMC::DMC(int n_p, int dim, int n_w, int n_c, int block_size, double E_T,
        Sampling *sampling,
        System *system,
        Kinetics *kinetics,
        Jastrow *jastrow,
        bool dist_from_file)
: QMC(n_p, dim, n_c, sampling, system, kinetics, jastrow) {

    this->dist_from_file = dist_from_file;

    this->block_size = block_size;
    this->n_w_orig = n_w;
    this->E_T = E_T;

}

void DMC::initialize() {

    E = 0;
    dmc_E = 0;

    jastrow->initialize();
    trial_walker = new Walker(n_p, dim);


    K = 10;
    int max_walkers = K * n_w_orig;

    Angry_mob = new Walker*[max_walkers];

    for (int i = 0; i < max_walkers; i++) {
        Angry_mob[i] = new Walker(n_p, dim);
    }


    if (dist_from_file) {

        ifstream dist;
        dist.open("dist_out.dat");

        for (int k = 0; k < n_w_orig; k++) {
            Angry_mob[k] = new Walker(n_p, dim);
            sampling->set_trial_pos(Angry_mob[k], true, &dist);

        }

        for (int k = n_w_orig; k < K * n_w_orig; k++) {
            Angry_mob[k] = new Walker(n_p, dim, false);
        }

        dist.close();

    } else {

        for (int k = 0; k < n_w_orig; k++) {

            Angry_mob[k] = new Walker(n_p, dim);

            while (Angry_mob[k]->is_singular()) {
                sampling->set_trial_pos(Angry_mob[k]);
            }

        }

        for (int k = n_w_orig; k < K * n_w_orig; k++) {
            Angry_mob[k] = new Walker(n_p, dim, false);
        }

    }

    n_w = n_w_orig;

}

void DMC::output() const {
    cout << "DMC energy: " << dmc_E / n_c << endl;
}

void DMC::Evolve_walker(double GB) {

    int branch_mean = int(GB + sampling->call_RNG()); //random int with mean=GB

    if (branch_mean == 0) {

        trial_walker->kill();

    } else {

        for (int n = 1; n < branch_mean; n++) {

            n_w += 1;
            copy_walker(trial_walker, Angry_mob[n_w - 1]);

        }

        E += GB * local_E;
        samples++;

    }
}

void DMC::bury_the_dead() {
    int index, index2;
    //cout << "called" << endl;
    int newborn = n_w - n_w_last;
    int j = 0;


    for (int i = 0; i < n_w; i++) {

        if (Angry_mob[i]->is_dead()) {

            /* Fills the newborn into the dead people's homes.
             * if j < newborn ensures that only newborn get dead people's homes
             * index = last born newborn
             * we then reset the newborn's array element.
             */

            if (j < newborn) {

                //Index of the last newborn
                index = newborn + n_w_last - 1 - j;

                copy_walker(Angry_mob[index], Angry_mob[i]);

                //                Angry_mob[index] = new Walker(n_p, dim, false); //THIS MEMORY EFFICIENT?

                j++;
            } else {

                index2 = n_w_last - 1;

                //Find the last alive walker
                while (Angry_mob[index2]->is_dead()) {
                    index2--;
                }

                /* test ensures previously moved is not moved again.
                 * compresses all alive to the start of the array.
                 */

                if (index2 > i) {
                    copy_walker(Angry_mob[index2], Angry_mob[i]);
                    Angry_mob[index2]->kill();
                }
            }
        }
    }

    //Finding the last alive walker below n_w_now
    for (int i = 0; i < n_w_last; i++) {
        if (!Angry_mob[i]->is_dead()) {
            ////
            ////            //clearing up if someone is dead;
            ////            Angry_mob[i] = new Walker(n_p, dim, false);
            //        } else {
            n_w = i + 1;
        }
    }

    //if there are more newborn than dead, at this point n_w = n_w_now;
    //adding the newborns which didn't get a home to that list:
    n_w += newborn - j;

}

void DMC::initialize_walker(int k) {
    original_walker = Angry_mob[k];
    copy_walker(original_walker, trial_walker);
    accepted = 0;
}

void DMC::update_energies() {
    dmc_E += E / samples;
    E_T = dmc_E / cycle; //+ (1. / 1000) * log((double) n_w_orig / n_w);
}

void DMC::run_method() {

    initialize();

    for (cycle = 1; cycle <= n_c; cycle++) {

        n_w_last = n_w;
        E = 0;
        samples = 0;

        for (int k = 0; k < n_w_last; k++) {

            initialize_walker(k);

            int b = 0;

            while ((!trial_walker->is_dead())*(b < block_size)) {

                calculate_energy_necessities(original_walker);
                double local_E_orig = calculate_local_energy(original_walker);

                diffuse_walker(original_walker, trial_walker);

                calculate_energy_necessities(trial_walker);
                local_E = calculate_local_energy(trial_walker);

                double GB = sampling->get_branching_Gfunc(local_E, local_E_orig, E_T);

                Evolve_walker(GB);

                b++;
            }
        }

        bury_the_dead();

        update_energies();


        ///DEBUGGING
        printf("%1.5f %1.5f %1.5f%%", dmc_E / cycle, (double) n_w / n_w_orig, (double) cycle / n_c * 100);
        cout << endl;

        dump_output();

    }

    finalize_output();
}

