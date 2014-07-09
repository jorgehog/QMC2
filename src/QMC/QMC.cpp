#include "QMC.h"

#include <iomanip>

#include "../structs.h"

#include "../Walker/Walker.h"
#include "../System/System.h"
#include "../Orbitals/Orbitals.h"
#include "../Sampling/Sampling.h"
#include "../Sampling/Importance/Importance.h"
#include "../Jastrow/Jastrow.h"
#include "../ErrorEstimator/ErrorEstimator.h"
#include "../OutputHandler/Distribution/Distribution.h"
#include "../Jastrow/No_Jastrow/No_Jastrow.h"


using namespace QMC2;

QMC::QMC(GeneralParams & gP, int n_c,
         SystemObjects & sO,
         ParParams & pp,
         bool silent,
         double dt,
         int n_w,
         int K) {

    n_p = gP.n_p;
    dim = gP.dim;
    this->n_c = n_c;
    this->n_w = n_w;
    n_w_size = K*n_w;

    dtOrig = dt;
    this->silent = silent;

    n2 = ceil(n_p / 2.0);

    trial_walker = new Walker(n_p, dim);
    original_walkers = new Walker*[n_w_size];
    for (unsigned int i = 0; i < n_w_size; i++) {
        original_walkers[i] = new Walker(n_p, dim);
    }

    if (sO.jastrow != NULL)
    {
        jastrow = sO.jastrow;
    }
    else
    {
        jastrow = new No_Jastrow();
    }

    if (sO.sample_method != NULL)
    {
        sampling = sO.sample_method;
    }
    else
    {
        sampling = new Importance(gP);
    }

    system = sO.system;

    sampling->set_qmc_ptr(this);
    get_orbitals_ptr()->set_qmc_ptr(this);

    accepted = 0;
    total_samples = 0;

    n_nodes = pp.n_nodes;
    node = pp.node;
    is_master = pp.is_master;
    parallel = pp.parallel;

    jastrow->initialize();

    runpath = gP.runpath;
    dist_path = runpath + "walker_positions/";
    last_inserted = 0;

    p_start = 0;
    if (gP.deadlock) {
        if (is_master) std::cout << "Initializing deadlock at " << gP.deadlock_x << std::endl;
        p_start = 1;
        sampling->set_deadlock(gP.deadlock_x); //Allowed since the constructor is a friend of Sampling
    }
}


void QMC::update_subsamples(double weight) {
    kinetic_sampler->update_mean(weight);
    system->update_potential_samples(weight);
}

void QMC::push_subsamples() {

    kinetic_sampler->push_mean();

    system->push_potential_samples();

    for (Sampler * sampler_method : samplers) {
        sampler_method->push_mean();
    }
}

void QMC::finalize()
{

    error_estimator->finalize();

    kinetic_sampler->getMeanErrorEstimator()->finalize();
    kinetic_sampler->getMeanOfMeansErrorEstimator()->finalize();

    system->finalize_potential_samples();

    for (Sampler * sampler_method : samplers) {
        sampler_method->getMeanErrorEstimator()->finalize();
        sampler_method->getMeanOfMeansErrorEstimator()->finalize();
    }
}


void QMC::reset_all()
{

    error_estimator->reset();

    kinetic_sampler->reset();

    system->reset_potential_samples();

    for (Sampler * sampler : samplers)
    {
        sampler->reset();
    }

    accepted = 0;
    total_samples = 0;
    last_inserted = 0;

    dist.reset();
    dist.set_size(dist_size, dim);

}

void QMC::initializeRun(const std::string method, int sampleType)
{
    sampling->set_dt(dtOrig);
    m_currentlyRunningMethod = method;

    kinetic_sampler->initializeErrorEstimator(sampleType, n_c);
    system->initializeErrorEstimator(sampleType, n_c);

    for (Sampler * sampler_method : samplers) {
        sampler_method->getMeanErrorEstimator()->setNumberOfCycles(n_c);
    }

}

void QMC::update_samplers(Walker *walker, double weight = 1.0)
{

    for (Sampler * sampler_method : samplers)
    {
        sampler_method->push_values(walker, weight);
    }

}


void QMC::dump_subsamples(bool mean_of_means) {

    using namespace std;

    stringstream s;

    s << "Kinetic " << setprecision(6) << fixed;
    if (mean_of_means) {
        s << kinetic_sampler->extract_mean_of_means();
    } else {
        s << kinetic_sampler->extract_mean();
    }

    s << endl;

    s << system->dump_samples(mean_of_means) << endl;

    if (is_master) {
        cout << setprecision(6) << fixed;
        cout << s.str();
    }

}

void QMC::finalize_distribution(){
    
    //scrap out all the over-allocated space (DMC)
    dist.resize(last_inserted, dim);

    if (dim == 3) {
        distribution->generate_distribution3D(dist, n_p);
    } else {
        distribution->generate_distribution2D(dist, n_p);
    }

    dist.reset();
}

void QMC::estimate_error()
{

    error = error_estimator->estimate_error();

    if (is_master) {
        if ((error != 0) && (!silent)) std::cout << "Estimated Error: " << error << std::endl;
    }

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
    double spatial_jast = walker_post->spatial_ratio * get_jastrow_ptr()->get_j_ratio(walker_post, walker_pre, particle);

    //    test_ratios(walker_pre, walker_post, particle, spatial_jast);

    double G = sampling->get_g_ratio(walker_post, walker_pre);

    return spatial_jast * spatial_jast * G;
}

void QMC::set_spin_state(unsigned int particle) const {
    unsigned int start = n2 * (particle >= n2);
    unsigned int end = start + n2 - 1;
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

    for (unsigned int i = 0; i < dim; i++) {
        walker_pre->r(particle, i) = walker_post->r(particle, i);
    }

    for (unsigned int i = 0; i < n_p; i++) {
        walker_pre->r_rel(i, particle) = walker_pre->r_rel(particle, i) = walker_post->r_rel(i, particle);

        for (unsigned int k = 0; k < dim; k++) {
            walker_pre->dJ(particle, i, k) = walker_post->dJ(particle, i, k);
            walker_pre->dJ(i, particle, k) = walker_post->dJ(i, particle, k);
        }
    }

    for (unsigned int i = 0; i < n2; i++) {
        walker_pre->phi(particle, i) = walker_post->phi(particle, i);
        walker_pre->dell_phi(particle)(arma::span(), i) = walker_post->dell_phi(particle)(arma::span(), i);
    }

    walker_pre->r2(particle) = walker_post->r2(particle);

    system->update_walker(walker_pre, walker_post, particle);
    sampling->update_walker(walker_pre, walker_post, particle);
}

void QMC::reset_walker(const Walker* walker_pre, Walker* walker_post, int particle) const {
    for (unsigned int i = 0; i < dim; i++) {
        walker_post->r(particle, i) = walker_pre->r(particle, i);
    }

    for (unsigned int i = 0; i < n_p; i++) {
        walker_post->r_rel(i, particle) = walker_post->r_rel(particle, i) = walker_pre->r_rel(i, particle);

        for (unsigned int k = 0; k < dim; k++) {
            walker_post->dJ(particle, i, k) = walker_pre->dJ(particle, i, k);
            walker_post->dJ(i, particle, k) = walker_pre->dJ(i, particle, k);
        }
    }

    for (unsigned int i = 0; i < n2; i++) {
        walker_post->phi(particle, i) = walker_pre->phi(particle, i);
        walker_post->dell_phi(particle)(arma::span(), i) = walker_pre->dell_phi(particle)(arma::span(), i);
    }

    walker_post->r2(particle) = walker_pre->r2(particle);

    system->reset_walker(walker_pre, walker_post, particle);
    sampling->reset_walker(walker_pre, walker_post, particle);
}



void QMC::diffuse_walker(Walker* original, Walker* trial) {
    for (unsigned int particle = p_start; particle < n_p; particle++) {

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
    unsigned int i, j;
    double xterm, e_kinetic;

    e_kinetic = xterm = 0;

    for (i = 0; i < n_p; i++) {
        for (j = 0; j < dim; j++) {
            xterm += walker->jast_grad(i, j) * walker->spatial_grad(i, j);
        }
    }

    e_kinetic = -xterm - 0.5 * walker->lapl_sum;

    kinetic_sampler->queue_value(e_kinetic);

    return e_kinetic;
}

void QMC::queue_weight(double weight)
{
    kinetic_sampler->queue_weight(weight);

    for (Potential *potential : system->potentials)
    {
        potential->
    }

}

void QMC::get_QF(Walker* walker) const {
    walker->qforce = 2 * (walker->jast_grad + walker->spatial_grad);
}

void QMC::get_accepted_ratio() {

    using namespace std;

    stringstream s;

    if (silent) return;

    double ratio = ((double) accepted) / total_samples;

    s << "Acceptance ratio: " << error_estimator->combine_mean(ratio, total_samples) << endl;

    if (is_master)
    {
        cout << s.str();
    }

}

void QMC::clean() {
    m_currentlyRunningMethod = "None";
    sampling->clear_deadlock();
}

/*

 Test functions
 *
 */

void QMC::test_gradients(Walker* walker) {

    double val = get_wf_value(walker);
    double h = 0.000001;
    for (unsigned int i = 0; i < n_p; i++) {
        for (unsigned int j = 0; j < dim; j++) {
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
    if ((abs(t1 - 1) < eps) & (abs(t2 - 1) < eps) & (abs(t3 - 1) < eps)) {
        std::cout << "ratio success" << std::endl;
    } else {
        std::cout << "ratio fail " << t1 << " " << t2 << " " << t3 << std::endl;
    }

}

//ADD

Orbitals* QMC::get_orbitals_ptr() const {
    return system->get_orbital_ptr();
}

double QMC::calculate_local_energy(const Walker* walker) {
    return get_KE(walker) + system->get_potential_energy(walker);
}

double QMC::get_wf_value(const Walker* walker) const {
    return system->get_spatial_wf(walker) * jastrow->get_val(walker);
}

void QMC::get_laplsum(Walker* walker) const {
    walker->lapl_sum = system->get_spatial_lapl_sum(walker) + jastrow->get_lapl_sum(walker);
}

std::string QMC::m_currentlyRunningMethod = "None";
