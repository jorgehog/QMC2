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

    double vmc_E, E2;

    Walker* original_walker;
    Walker* trial_walker;

    virtual void initialize();

    virtual bool move_autherized(double A);

    void calculate_energy(Walker* walker);
    void scale_values();

public:

    VMC(GeneralParams &, VMCparams &, SystemObjects &);

    double get_var() const;
    double get_energy() const;
    double get_e2() const;
    void set_e(double e);
    void set_e2(double e2);

    virtual void run_method();
    virtual void user_output() const;

    friend class Minimizer;
    friend class ASGD;

    friend class Distribution;
    friend class BlockingData;

};

#endif	/* VMC_H */
