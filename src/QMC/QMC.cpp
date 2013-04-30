/* 
 * File:   QMC.cpp
 * Author: jorgehog
 * 
 * Created on 30. mars 2012, 17:42
 */

#include "../QMCheaders.h"

QMC::QMC(GeneralParams & gP, int n_c,
        SystemObjects & sO,
        ParParams & pp,
        int n_w,
        int K) {

    n_p = gP.n_p;
    dim = gP.dim;
    this->n_c = n_c;
    this->n_w = n_w;
    n_w_size = K*n_w;

    n2 = n_p / 2;

    trial_walker = new Walker(n_p, dim);
    original_walkers = new Walker*[n_w_size];
    for (int i = 0; i < n_w_size; i++) {
        original_walkers[i] = new Walker(n_p, dim);
    }

    jastrow = sO.jastrow;
    sampling = sO.sample_method;
    system = sO.SYSTEM;

    sampling->set_qmc_ptr(this);
    get_orbitals_ptr()->set_qmc_ptr(this);

    accepted = 0;
    total_samples = 0;

    n_nodes = pp.n_nodes;
    node = pp.node;
    is_master = pp.is_master;
    parallel = pp.parallel;

    jastrow->initialize();

    if (is_master) {
        std_out = new STDOUT();
    } else {
        std_out = new NO_STDOUT();
    }

    runpath = gP.runpath;
    dist_path = runpath + "walker_positions/";
    last_inserted = 0;

}

QMC::QMC() {

}

void QMC::update_subsamples(double weight) {
    kinetic_sampler.update_mean(weight);
    system->update_potential_samples(weight);
}

void QMC::push_subsamples(){
    kinetic_sampler.push_mean();
    system->push_potential_samples();
}

void QMC::dump_subsamples(bool mean_of_means) {

    using namespace std;

    s << "Kinetic " << setprecision(6) << fixed;
    if (mean_of_means){
        s << kinetic_sampler.extract_mean_of_means();
    } else {
        s << kinetic_sampler.extract_mean();
    }
    
    s << endl;
    
    s << system->dump_samples(mean_of_means) << endl;

    if (is_master) {
        cout << setprecision(6) << fixed;
        cout << s.str();
    }

    s.str(string());

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

void QMC::estimate_error() const {

    double error = error_estimator->estimate_error();

    if (is_master) {
        if (error != 0) std::cout << "Estimated Error: " << error << std::endl;
    }

    error_estimator->finalize();

}

void QMC::get_gradients(const Walker* walker_pre, Walker* walker_post, int particle) const {
    jastrow->get_grad(walker_pre, walker_post, particle);
    system->get_spatial_grad(walker_post, particle);
}

void QMC::get_gradients(Walker* walker) const {
    jastrow->get_grad(walker);
    system->get_spatial_grad_full(walker);
}

double QMC::get_acceptance_ratio(const Walker* walker_pre, const Walker* walker_post, int particle) const {
    double spatial_jast = sampling->get_spatialjast_ratio(walker_post, walker_pre, particle);

    //    test_ratios(walker_pre, walker_post, particle, spatial_jast);

    double G = sampling->get_g_ratio(walker_post, walker_pre);

    return spatial_jast * spatial_jast * G;
}

void QMC::set_spin_state(int particle) const {
    int start = n2 * (particle >= n2);
    int end = start + n2 - 1;
    system->set_spin_state(start, end);
    sampling->set_spin_state(start, end);
}

bool QMC::metropolis_test(double A) {
    double r = sampling->call_RNG();

    if (r <= A) {
        //        if (!system->allow_transition()) {
        //            std::cout << "move " << total_samples << " rejected due to fixed node" << std::endl;
        //        }
        return true;

    } else {
        return false;
    }
}

void QMC::calculate_energy_necessities(Walker* walker) const {
    sampling->calculate_energy_necessities(walker);
    get_laplsum(walker);
}

void QMC::update_walker(Walker* walker_pre, const Walker* walker_post, int particle) const {

    for (int i = 0; i < dim; i++) {
        walker_pre->r(particle, i) = walker_post->r(particle, i);
    }

    for (int i = 0; i < n_p; i++) {
        walker_pre->r_rel(i, particle) = walker_pre->r_rel(particle, i) = walker_post->r_rel(i, particle);

        for (int k = 0; k < dim; k++) {
            walker_pre->dJ(particle, i, k) = walker_post->dJ(particle, i, k);
            walker_pre->dJ(i, particle, k) = walker_post->dJ(i, particle, k);
        }
    }

    for (int i = 0; i < n2; i++) {
        walker_pre->phi(particle, i) = walker_post->phi(particle, i);
        walker_pre->dell_phi(particle)(arma::span(), i) = walker_post->dell_phi(particle)(arma::span(), i);
    }

    walker_pre->r2(particle) = walker_post->r2(particle);

    system->update_walker(walker_pre, walker_post, particle);
    sampling->update_walker(walker_pre, walker_post, particle);
}

void QMC::reset_walker(const Walker* walker_pre, Walker* walker_post, int particle) const {
    for (int i = 0; i < dim; i++) {
        walker_post->r(particle, i) = walker_pre->r(particle, i);
    }

    for (int i = 0; i < n_p; i++) {
        walker_post->r_rel(i, particle) = walker_post->r_rel(particle, i) = walker_pre->r_rel(i, particle);

        for (int k = 0; k < dim; k++) {
            walker_post->dJ(particle, i, k) = walker_pre->dJ(particle, i, k);
            walker_post->dJ(i, particle, k) = walker_pre->dJ(i, particle, k);
        }
    }

    for (int i = 0; i < n2; i++) {
        walker_post->phi(particle, i) = walker_pre->phi(particle, i);
        walker_post->dell_phi(particle)(arma::span(), i) = walker_pre->dell_phi(particle)(arma::span(), i);
    }

    walker_post->r2(particle) = walker_pre->r2(particle);

    system->reset_walker(walker_pre, walker_post, particle);
    sampling->reset_walker(walker_pre, walker_post, particle);
}

void QMC::diffuse_walker(Walker* original, Walker* trial) {
    for (int particle = 0; particle < n_p; particle++) {

        set_spin_state(particle);
        sampling->update_pos(original, trial, particle);

        //        test_gradients(trial);

        double A = get_acceptance_ratio(original, trial, particle);

        if (move_autherized(A)) {
            accepted++;
            update_walker(original, trial, particle);
        } else {
            reset_walker(original, trial, particle);
        }
        total_samples++;

    }
}

void QMC::copy_walker(const Walker* parent, Walker* child) const {

    child->r2 = parent->r2;
    child->r = parent->r;
    child->r_rel = parent->r_rel;

    child->phi = parent->phi;
    child->dell_phi = parent->dell_phi;
    child->dJ = parent->dJ;

    child->ressurect();
    child->set_E(parent->get_E());
    system->copy_walker(parent, child);
    sampling->copy_walker(parent, child);

}

double QMC::get_KE(const Walker* walker) {
    int i, j;
    double xterm, e_kinetic;

    e_kinetic = xterm = 0;

    for (i = 0; i < n_p; i++) {
        for (j = 0; j < dim; j++) {
            xterm += walker->jast_grad(i, j) * walker->spatial_grad(i, j);
        }
    }

    e_kinetic = -xterm - 0.5 * walker->lapl_sum;

    kinetic_sampler.queue_value(e_kinetic);

    return e_kinetic;
}

void QMC::get_QF(Walker* walker) const {
    walker->qforce = 2 * (walker->jast_grad + walker->spatial_grad);
}

void QMC::get_accepted_ratio() {

    double ratio = ((double) accepted) / total_samples;

    s << "Acceptance ratio: " << error_estimator->combine_mean(ratio, total_samples);
    s << std::endl;
    std_out->cout(s);
}

/*
 
 Test functions
 * 
 */

void QMC::test_gradients(Walker* walker) {

    double val = get_wf_value(walker);
    double h = 0.000001;
    for (int i = 0; i < n_p; i++) {
        for (int j = 0; j < dim; j++) {
            walker->r(i, j) += h;
            walker->make_rel_matrix();
            walker->calc_r_i2();
            sampling->set_trial_states(walker);

            double vp = get_wf_value(walker);

            walker->r(i, j) -= 2 * h;
            walker->make_rel_matrix();
            walker->calc_r_i2();
            sampling->set_trial_states(walker);

            double vm = get_wf_value(walker);

            double deriv = (vp - vm) / (2 * val * h);
            double anal = walker->spatial_grad(i, j) + walker->jast_grad(i, j);
            if ((deriv - anal) / anal > 1E-3) {
                std::cout << i << "  " << j << "   " << deriv << "   " << anal << std::endl;
            }
            walker->r(i, j) += h;
            walker->make_rel_matrix();
            walker->calc_r_i2();
            sampling->set_trial_states(walker);
        }
    }


}

void QMC::test_ratios(const Walker* walker_pre, const Walker* walker_post, int particle, double R_qmc) const {

    double wf_new_s = system->get_spatial_wf(walker_post);
    double wf_new_j = jastrow->get_val(walker_post);
    double wf_new = wf_new_s*wf_new_j;

    double wf_old_s = system->get_spatial_wf(walker_pre);
    double wf_old_j = jastrow->get_val(walker_pre);
    double wf_old = wf_old_s*wf_old_j;


    double R_opt_s = system->get_spatial_ratio(walker_pre, walker_post, particle);
    double R_opt_j = jastrow->get_j_ratio(walker_post, walker_pre, particle);
    double R = wf_new / wf_old;

    double t1 = R_opt_s / (wf_new_s / wf_old_s);
    double t2 = R_opt_j / (wf_new_j / wf_old_j);
    double t3 = R_qmc / R;

    double eps = 1e-10;
    if (abs(t1 - 1) < eps & abs(t2 - 1) < eps & abs(t3 - 1) < eps) {
        std::cout << "ratio success" << std::endl;
    } else {
        std::cout << "ratio fail " << t1 << " " << t2 << " " << t3 << std::endl;
    }

}