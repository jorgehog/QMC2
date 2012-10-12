/* 
 * File:   Jastrow.h
 * Author: jorgehog
 *
 * Created on 30. mars 2012, 16:52
 */

#ifndef JASTROW_H
#define	JASTROW_H

class Jastrow {
protected:
    int n_p;
    int n2;
    int dim;
    
    bool active;

    virtual double get_parameter(int n) = 0;
    virtual void set_parameter(double param, int n) = 0;
    virtual double get_variational_derivative(const Walker* walker, int n) const = 0;


public:
    Jastrow(int n_p, int dim);
    Jastrow();

    virtual void initialize() = 0;

    virtual double get_val(const Walker* walker) const = 0;
    virtual double get_j_ratio(const Walker* walker_new, const Walker* walker_old, int i) const = 0;
    virtual void get_grad(Walker* walker) const = 0;
    virtual double get_lapl_sum(const Walker* walker) const = 0;
  
    friend class Minimizer;
    friend class ASGD;
    
};

#endif	/* JASTROW_H */
