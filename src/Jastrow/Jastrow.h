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
    virtual double get_variational_derivative(const Walker* walker, int n);
    
    double get_derivative_num(Walker* walker, int i, int d) const;
    double get_laplaciansum_num(Walker* walker) const;

public:
    Jastrow(int n_p, int dim);
    Jastrow();

    virtual void initialize() = 0;

    virtual double get_val(const Walker* walker) const = 0;
    virtual double get_j_ratio(const Walker* walker_new, const Walker* walker_old, int i) const = 0;
    
    virtual void get_grad(Walker* walker) const = 0;
    virtual void get_grad(const Walker* walker_pre, Walker* walker_post, int i) const = 0;
    virtual void get_dJ_matrix(Walker* walker, int i) const = 0;
    void get_dJ_matrix(Walker* walker) const;
    
    virtual double get_lapl_sum(Walker* walker) const = 0;
    
  
    friend class Minimizer;
    friend class ASGD;
    friend class stdoutASGD;
    
};

#endif	/* JASTROW_H */
