/* 
 * File:   QMC.cpp
 * Author: jorgehog
 * 
 * Created on 30. mars 2012, 17:42
 */

#include "../QMCheaders.h"

QMC::QMC(int n_p, int dim, int n_c, Jastrow *jastrow, Sampling *sampling, System *system, Kinetics *kinetics) {

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

    update_necessities(original_walker, trial_walker, particle);

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

void QMC::diffuse_walker() {
    for (int particle = 0; particle < n_p; particle++) {

        sampling->reset_control_parameters();
        update_pos(original_walker, trial_walker, particle);

        double A = get_acceptance_ratio(original_walker, trial_walker, particle);

        if (metropolis_test(A)) {
            update_walker(original_walker, trial_walker, particle);
        } else {
            reset_walker(original_walker, trial_walker, particle);
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

VMC::VMC(int n_p, int dim, int n_c, Jastrow *jastrow, Sampling *sampling, System *system, Kinetics *kinetics, bool dist_to_file)
: QMC(n_p, dim, n_c, jastrow, sampling, system, kinetics) {

    vmc_E = 0;
    E2 = 0;

    original_walker = new Walker(n_p, dim);
    trial_walker = new Walker(n_p, dim);

    this->dist_to_file = dist_to_file;
}

void VMC::initialize() {
    jastrow->initialize();

    sampling->reset_control_parameters();
    sampling->set_trial_pos(original_walker);

    copy_walker(original_walker, trial_walker);
}

void VMC::calculate_energy(Walker* walker) {

    double local_E = calculate_local_energy(walker);

    vmc_E += local_E;
    E2 += local_E*local_E;
}

void VMC::scale_values() {
    vmc_E /= n_c;
    E2 /= n_c;
}

void VMC::run_method() {

    //WRAP INTO OBJECT CLASS OUTPUT.init
    ofstream dist;
    dist.open("dist_out.dat");
    //

    initialize();

    //OK
    //    cout << original_walker->r << original_walker->r_rel << original_walker->r2 << trial_walker->value << endl;
    //    cout << "-----------------------------------------\n";
    //    cout << trial_walker->r << trial_walker->r_rel << trial_walker->r2 << trial_walker->value << endl;
    //    exit(1);
    double ek = 0;
    double ep = 0;
    for (int cycle = 1; cycle <= n_c + thermalization; cycle++) {

        diffuse_walker();

        //WRAP INTO OBJECT CLASS OUTPUT.output
        //        if ((cycle > n_c / 2) && (cycle % 100 == 0)) {
        //            for (int i = 0; i < n_p; i++) {
        //                for (int j = 0; j < dim; j++) {
        //                    if (j == dim - 1) {
        //                        dist << original_walker->r(i, j);
        //                    } else {
        //                        dist << original_walker->r(i, j) << " ";
        //                    }
        //                }
        //                dist << endl;
        //            }
        //        }
        ///


        if (cycle > thermalization) {
            ek += sum(original_walker->r2);
            calculate_energy_necessities(original_walker);
            calculate_energy(original_walker);
        }
    }


    ek /= n_c;
    ep /= n_c;
    cout << ek << "\t" << ep << endl;
    scale_values();

    //WRAP INTO OBJECT CLASS OUTPUT.finalize
    dist.close();

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

DMC::DMC(int n_p, int dim, int n_w, int n_c, double E_T, Jastrow *jastrow, Sampling *sampling,
        System *system, Kinetics *kinetics, bool dist_from_file)
: QMC(n_p, dim, n_c, jastrow, sampling, system, kinetics) {

    this->dist_from_file = dist_from_file;
    this->n_w_orig = n_w;
    this->E_T = E_T;
    E = 0;
}

void DMC::initialize() {

    jastrow->initialize();
    original_walker = new Walker(n_p, dim);

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

    E = E_T;
    n_w = n_w_orig;
}

//void DMC::increase_walker_space() {
//    int i, j, k, l;
//
//    //creating tmp-replacement
//    Walker **tmp = new Walker*[2 * K * n_w_orig];
//
//    for (i = 0; i < 2 * K * n_w_orig; i++) {
//        tmp[i] = new Walker(this);
//    }
//
//    //copying current walkers into the tmp
//    for (l = 0; l < K * n_w_orig; l++) {
//        for (i = 0; i < n_p; i++) {
//            for (j = 0; j < dim; j++) {
//                tmp[l]->r[i][j] = Angry_mob[l]->r[i][j];
//            }
//
//            for (k = i + 1; k < n_p; k++) {
//                tmp[l]->r_rel[i][k] = tmp[l]->r_rel[k][i] = Angry_mob[l]->r[i][k];
//            }
//        }
//
//        if (is_importance()) {
//            for (i = 0; i < n_p; i++) {
//                for (j = 0; j < dim; j++) {
//                    tmp[l]->qforce[i][j] = Angry_mob[l]->qforce[i][j];
//                }
//            }
//        }
//        tmp[l]->is_murdered = Angry_mob[l]->is_murdered;
//    }
//
//    //deleting the old
//    for (i = 0; i < K * n_w_orig; i++) {
//        delete Angry_mob[i];
//    }
//
//    delete[] Angry_mob;
//
//    //creating a new with more room
//    K *= 2;
//    Walker** Angry_mob = new Walker*[K * n_w_orig];
//
//    for (i = 0; i < K * n_w_orig; i++) {
//        Angry_mob[i] = new Walker(this);
//    }
//
//
//
//    //copying from tmp into the new
//    for (l = 0; l < K * n_w_orig; l++) {
//        for (i = 0; i < n_p; i++) {
//            for (j = 0; j < dim; j++) {
//                Angry_mob[l]->r[i][j] = tmp[l]->r[i][j];
//            }
//
//            for (k = i + 1; k < n_p; k++) {
//                Angry_mob[l]->r_rel[i][k] = Angry_mob[l]->r_rel[k][i] = tmp[l]->r[i][k];
//            }
//        }
//
//        if (is_importance()) {
//            for (i = 0; i < n_p; i++) {
//                for (j = 0; j < dim; j++) {
//                    Angry_mob[l]->qforce[i][j] = tmp[l]->qforce[i][j];
//                }
//            }
//        }
//        Angry_mob[l]->is_murdered = tmp[l]->is_murdered;
//
//    }
//
//
//    //deleting tmp
//    for (i = 0; i < K * n_w_orig; i++) {
//        delete tmp[i];
//    }
//
//    delete[] tmp;
//
//    this->Angry_mob = Angry_mob;
//
//}

void DMC::output() const {
    cout << "DMC energy: " << E / n_c << endl;
}

void DMC::Evolve_walker(double GB) {
    int n;
    int branch_mean;

    branch_mean = int(GB + sampling->call_RNG()); //random int with mean=GB

    if (branch_mean >= 1) {
        for (n = 0; n < branch_mean - 1; n++) {

            n_w += 1;
            copy_walker(trial_walker, Angry_mob[n_w - 1]);
        }

        //calculate_energy_necessities(trial_walker);
        //get_wf_value(trial_walker);
        //E_avg += GB * calculate_local_energy(trial_walker);
        //samples++;

    } else if (branch_mean == 0) {
        trial_walker->kill();
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

void DMC::increase_walker_space() {

}

void DMC::initialize_walker(int k) {
    trial_walker = Angry_mob[k];
    copy_walker(trial_walker, original_walker);
    accepted = 0;
}

void DMC::update_energies(int n) {
    //E += E_avg / samples;
    E_T = E / (n + 1) + (1. / 1000) * log((double) n_w_orig / n_w);
}

void DMC::run_method() {

    initialize();

    //DEBUGGING
    ofstream file;
    file.open("e_vs_t.dat");


    for (int cycle = 0; cycle < n_c; cycle++) {
        n_w_last = n_w;
        E_avg = 0;
        samples = 0;
        for (int k = 0; k < n_w_last; k++) {
            initialize_walker(k);


            for (int particle = 0; particle < n_p; particle++) {
                update_pos(original_walker, trial_walker, particle);
                //metro+++
            }

            double GB = sampling->get_branching_Gfunc(original_walker, trial_walker, E_T);
            cout << GB << endl;
            Evolve_walker(GB);
        }

        bury_the_dead();

        update_energies(cycle);


        ///DEBUGGING
        if (cycle % (n_c / 1000) == 0) {
            cout << "\r";
            printf("%1.5f %1.5f %1.5f %1.5f%%", E_T, E / (cycle + 1), (double) n_w / n_w_orig, (double) (cycle + 1) / n_c * 100);
            file << (cycle + 1) * 0.01 << "\t" << E / (cycle + 1) << "\t" << (double) n_w / n_w_orig << endl;
        }

    }
}

