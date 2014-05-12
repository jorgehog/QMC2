#include "Sampling.h"

#include "../Walker/Walker.h"
#include "../System/System.h"
#include "../Jastrow/Jastrow.h"
#include "../Orbitals/Orbitals.h"


using namespace QMC2;

Sampling::Sampling(int n_p, int dim) {

    this->n_p = n_p;
    this->dim = dim;

    n2 = ceil(n_p / 2.0);
    deadlock = false;
}

Sampling::Sampling() {

}

void Sampling::set_deadlock(const double deadlock_x){
    
    if (deadlock) {
        std::cout << "Overriding deadlock detected. Exiting..." << std::endl;
        exit(1);
    }
    
    deadlock = true;
    this->deadlock_x = deadlock_x;
}

void Sampling::clear_deadlock() {
    deadlock = false;
}

void Sampling::set_trial_pos(Walker* walker) {

    for (int i = 0; i < n_p; i++) {
        for (int j = 0; j < dim; j++) {
            walker->r(i, j) = diffusion->Diffusion::get_new_pos(walker, i, j);
        }
    }

    //If test on deadlocks. Inefficient, however, trial positions are only set once in VMC.
    if (deadlock) {
        walker->r(0, 0) = deadlock_x;
        
        for (int j = 1; j < dim; j++){
            walker->r(0, j) = 0;
        }
    }
    
    walker->calc_r_i2();
    
    set_trial_states(walker);

    qmc->get_system_ptr()->initialize(walker);

    walker->make_rel_matrix();

    qmc->get_jastrow_ptr()->get_dJ_matrix(walker);

    get_necessities(walker);

}

void Sampling::set_trial_states(Walker * walker) {

    for (int i = 0; i < n_p; i++) {

        qmc->get_orbitals_ptr()->set_qnum_indie_terms(walker, i);

        for (int j = 0; j < n2; j++) {
            walker->phi(i, j) = qmc->get_orbitals_ptr()->phi(walker, i, j);
        }

        for (int j = 0; j < n2; j++) {
            for (int k = 0; k < dim; k++) {
                walker->dell_phi(i)(k, j) = qmc->get_orbitals_ptr()->del_phi(walker, i, j, k);
            }
        }
    }
}

void Sampling::update_pos(const Walker* walker_pre, Walker* walker_post, int particle) const {

    for (int j = 0; j < dim; j++) {
        walker_post->r(particle, j) = walker_pre->r(particle, j)
                + diffusion->get_new_pos(walker_pre, particle, j);
    }

    for (int j = 0; j < n_p; j++) {
        if (j != particle) {
            walker_post->r_rel(particle, j) = walker_post->r_rel(j, particle)
                    = walker_post->calc_r_rel(particle, j);
        }
    }

    walker_post->calc_r_i2(particle);
    qmc->get_orbitals_ptr()->set_qnum_indie_terms(walker_post, particle);

    for (int j = 0; j < n2; j++) {
        walker_post->phi(particle, j) = qmc->get_orbitals_ptr()->phi(walker_post, particle, j);
    }

    for (int j = 0; j < n2; j++) {
        for (int k = 0; k < dim; k++) {
            walker_post->dell_phi(particle)(k, j) = qmc->get_orbitals_ptr()->del_phi(walker_post, particle, j, k);
        }
    }

    qmc->get_jastrow_ptr()->get_dJ_matrix(walker_post, particle);
    walker_post->spatial_ratio = qmc->get_system_ptr()->get_spatial_ratio(walker_pre, walker_post, particle);
    qmc->get_system_ptr()->calc_for_newpos(walker_pre, walker_post, particle);

    update_necessities(walker_pre, walker_post, particle);

}
