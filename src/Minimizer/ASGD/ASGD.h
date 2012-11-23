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
    
    int sample;

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

    arma::rowvec parameter;

    arma::rowvec gradient;
    arma::rowvec gradient_local;

    arma::rowvec gradient_old;
    arma::rowvec gradient_tot;

    double f(double x) {
        return f_min + (f_max - f_min) / (1 - (f_max / f_min) * exp(-x / w));
    }
    
    void get_variational_gradients(Walker* walker, double e_local);

public:
    ASGD(VMC*, MinimizerParams &, const ParParams &);

    virtual VMC* minimize();

    friend class stdoutASGD;
};


#endif	/* ASGD_H */