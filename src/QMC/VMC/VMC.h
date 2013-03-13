/* 
 * File:   VMC.h
 * Author: jorgmeister
 *
 * Created on October 12, 2012, 2:43 PM
 */

#ifndef VMC_H
#define	VMC_H

class VMC : public QMC {
protected:

    int pop_tresh;
    int offset;
    int last_walker;
    
    double vmc_E;
    Walker* original_walker;

    void set_trial_positions();
    
    void store_walkers();
    void save_distribution();
    
    bool move_autherized(double A) {
        return metropolis_test(A);
    }

    void scale_values() {
        vmc_E /= (n_c*n_nodes);
    }
    
    void node_comm();

public:

    VMC(GeneralParams &, VMCparams &, SystemObjects &, ParParams &, int n_w, bool dist_out);
    ~VMC();
    
    void set_e(double E) {
        vmc_E = E;
    }

    double get_energy() const {
        return vmc_E;
    }

    void run_method();
    void output();

    friend class DMC;
    
    friend class Minimizer;
    friend class ASGD;

    friend class BlockingData;

};

#endif	/* VMC_H */
