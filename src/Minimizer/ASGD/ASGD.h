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

    rowvec parameter;
    
    rowvec gradient;
    rowvec gradient_local;
    
    rowvec gradient_old;
    rowvec gradient_tot;

    double f(double x);
    void get_variational_gradients(double e_local);

public:
    ASGD(VMC* vmc, int SDSsamples, int n_walkers, int n_c, int n_c_SGD, double f_min, double f_max, double w,
            double a = 1, double A = 1, int NSP = 1, int NJP = 1);

    virtual VMC* minimize();

};


#endif	/* ASGD_H */