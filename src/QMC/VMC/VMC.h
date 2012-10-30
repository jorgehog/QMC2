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

    double vmc_E;

    Walker* original_walker;
    Walker* trial_walker;

    virtual void initialize();

    virtual bool move_autherized(double A) {
        return metropolis_test(A);
    };

    void scale_values() {
        vmc_E /= n_c;
    };

public:

    VMC(GeneralParams &, VMCparams &, SystemObjects &);

    void set_e(double E) {
        vmc_E = E;
    }

    double get_energy() const {
        return vmc_E;
    }

    virtual void run_method();
    virtual void user_output() const;

    friend class Minimizer;
    friend class ASGD;

    friend class Distribution;
    friend class BlockingData;

};

#endif	/* VMC_H */
