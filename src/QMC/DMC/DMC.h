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

    int n_w;
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

    bool dist_from_file;
    std::string dist_in_path;

    Walker **original_walkers;
    Walker *trial_walker;

    void initialize();

    virtual bool move_autherized(double A) {
        return metropolis_test(A)&(A > 0);
    };

    void iterate_walker(int k, int n_b = 1, bool production = false);

    void Evolve_walker(int k, double GB);

    void bury_the_dead();

    void update_energies();

    void reset_parameters() {
        n_w_last = n_w;
        E = samples = deaths = 0;
    }

    virtual void node_comm();

    void switch_souls(int root, int root_id, int dest, int dest_id);

    void normalize_population();

public:

    static const int K = 2; //Factor of empty space for walkers over initial walkers
    static const int check_thresh = 5;
    static const int sendcount_thresh = 20;

    DMC(GeneralParams &, DMCparams &, SystemObjects &, ParParams &);

    virtual void run_method();

    virtual void output();

    friend class stdoutDMC;

};

#endif	/* DMC_H */
