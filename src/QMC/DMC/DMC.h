/* 
 * File:   DMC.h
 * Author: jorgmeister
 *
 * Created on October 12, 2012, 2:42 PM
 */

#ifndef DMC_H
#define	DMC_H

class DMC : public QMC {
protected:

    int n_w_last;
    
    int n_w_tot;
    arma::uvec n_w_list;
    
    int deaths;

    int block_size;
    int samples;
    unsigned long tot_samples;
    double E_tot;

    double dmc_E;
    double dmc_E_unscaled;
    double E_T;
    double E;

    void set_trial_positions();

    bool move_autherized(double A) {
        return metropolis_test(A)&system->allow_transition();
    };

    void iterate_walker(int k, int n_b = 1);

    void Evolve_walker(int k, double GB);

    void bury_the_dead();

    void update_energies();

    void reset_parameters() {
        n_w_last = n_w;
        E = samples = deaths = 0;
    }
    
    void node_comm();
    
    void save_distribution();

    void switch_souls(int root, int root_id, int dest, int dest_id);

    void normalize_population();
    

public:

    static const int K = 50; //Factor of empty space for walkers over initial walkers
    static const int check_thresh = 25; //normalize every [] cycle
    static const int sendcount_thresh = 20; //abort normalization if max difference in walkers small

    DMC(GeneralParams &, DMCparams &, SystemObjects &, ParParams &, VMC* vmc, bool dist_out);

    void run_method();

    void output();

    friend class stdoutDMC;

};

#endif	/* DMC_H */
