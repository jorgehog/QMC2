/* 
 * File:   QMC.h
 * Author: jorgehog
 *
 * Created on 30. mars 2012, 17:42
 */

#ifndef QMC_H
#define	QMC_H

#include "../QMCheaders.h"

/*! \brief The QMC superduper class!
 * 
 * This class is so cool!
 * 
 */
class QMC {
protected:

    STDOUT* std_out;
    std::stringstream s; //!< This stream is awesome!
    std::string name;
    std::string runpath;
    std::string dist_path;
    
    static const int n_distout = 100000;

    bool is_master;
    bool parallel;
    int node;
    int n_nodes;

    int dist_tresh;

    int n_c;
    int n_w;

    int n_p;
    int n2;
    int dim;

    int cycle;

    int accepted;
    int total_samples;
    int thermalization;

    double local_E;

    Walker *trial_walker;
    Walker **original_walkers;

    Jastrow *jastrow;
    Sampling *sampling;
    System *system;
    ErrorEstimator* error_estimator;

    std::vector<OutputHandler*> output_handler;

    virtual void set_trial_positions() = 0;

    virtual bool move_autherized(double A) = 0;

    void dump_output();
    void finalize_output();
    void dump_distribution(Walker* walker, std::string suffix);

    void diffuse_walker(Walker* original, Walker* trial);

    double get_acceptance_ratio(const Walker* walker_pre, const Walker* walker_post, int particle) const;

    void set_spin_state(int particle) const;

    bool metropolis_test(double A);
    void update_walker(Walker*walker_pre, const Walker* walker_post, int particle) const;
    void reset_walker(const Walker* walker_pre, Walker* walker_post, int particle) const;

    void copy_walker(const Walker* parent, Walker* child) const;

    void calculate_energy_necessities(Walker* walker) const;

    void estimate_error() const;

    virtual void node_comm() = 0;
    
    virtual void store_walkers() = 0;

    void switch_souls(int root, int source);

    //TEST FUNCTIONS
    void test_ratios(const Walker* walker_pre, const Walker* walker_post, int particle, double R_qmc) const;

public:

    QMC(GeneralParams &, int n_c,
            SystemObjects &,
            ParParams &,
            int K = 1);
    QMC();

    void add_output(OutputHandler* output_handler);

    virtual void run_method() = 0;
    virtual void output() = 0;



    double get_KE(const Walker* walker) const;
    void get_QF(Walker* walker) const;

    void get_gradients(const Walker* walker_pre, Walker* walker_post, int particle) const;
    void get_gradients(Walker* walker) const;

    void get_laplsum(Walker* walker) const {
        walker->lapl_sum = system->get_spatial_lapl_sum(walker) + jastrow->get_lapl_sum(walker);
    }

    double get_wf_value(const Walker* walker) const {
        return system->get_spatial_wf(walker) * jastrow->get_val(walker);
    }

    double calculate_local_energy(const Walker* walker) const;

    System* get_system_ptr() const {
        return system;
    }

    Sampling* get_sampling_ptr() const {
        return sampling;
    }

    Jastrow* get_jastrow_ptr() const {
        return jastrow;
    }

    Orbitals* get_orbitals_ptr() const {
        return system->get_orbital_ptr();
    }

    void get_accepted_ratio();

    void set_error_estimator(ErrorEstimator* error_estimator) {
        this->error_estimator = error_estimator;
    }

};

#endif	/* QMC_H */
