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
    ASGD(VMC*, MinimizerParams &);

    virtual VMC* minimize();
    virtual VMC* minimizeTEST();
    double TESTWF(Walker* walker);
    double TEST_E(Walker* walker);
    double TEST_G(Walker* walker_post, Walker* walker_pre);
    void TEST_DIFF(Walker* original, Walker* trial);
    std::ofstream DEBAG;
    long random_seed;

};


#endif	/* ASGD_H */