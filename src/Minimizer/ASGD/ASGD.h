/* 
 * File:   ASGD.h
 * Author: jorgehog
 *
 * Created on 23. August 2012, 16:52
 */

#ifndef ASGD_H
#define	ASGD_H

class ASGD : public Minimizer {
protected:

    int n_c;
    int n_c_SGD;
    int SGDsamples;
    int n_walkers;
    int thermalization;

    int k;

    double t_prev;
    double t;
    double step;
    double max_step;

    double E;

    double a;
    double A;
    double f_min;
    double f_max;
    double w;

    Walker** walkers;
    Walker** trial_walkers;

    rowvec parameter;

    rowvec gradient;
    rowvec gradient_local;

    rowvec gradient_old;
    rowvec gradient_tot;

    double f(double x);
    void get_variational_gradients(Walker* walker, double e_local);

public:
    ASGD(VMC* vmc,
            const rowvec & alpha,
            const rowvec & beta,
            int SGDsamples,
            int n_walkers,
            int n_c,
            int thermalization,
            int n_c_SGD,
            double max_step,
            double f_min,
            double f_max,
            double w,
            double a = 1,
            double A = 1);

    virtual VMC* minimize();
    std::ofstream DEBAG;

};


#endif	/* ASGD_H */