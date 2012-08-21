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


public:
    Jastrow(int n_p, int dim);
    Jastrow();

    virtual void initialize() = 0;

    virtual double get_val(const Walker* walker) const = 0;
    virtual double get_j_ratio(const Walker* walker_new, const Walker* walker_old, int i) const = 0;
    virtual void get_grad(Walker* walker) const = 0;
    virtual double get_lapl_sum(const Walker* walker) const = 0;

};

class No_Jastrow : public Jastrow {
public:

    No_Jastrow();

    virtual void get_grad(Walker* walker) const{

    };

    virtual void initialize() {
    };

    virtual double get_j_ratio(const Walker* walker_post, const Walker* walker_pre, int i) const {
        return 1;
    };

    virtual double get_val(const Walker* walker) const {
        return 1;
    };

    virtual double get_lapl_sum(const Walker* walker) const {
        return 0;
    };
};

class Pade_Jastrow : public Jastrow {
protected:
    double beta;
    arma::mat a;

public:

    Pade_Jastrow(int n_p, int dim, double beta);

    virtual void initialize();

    virtual void get_grad(Walker* walker) const;

    virtual double get_j_ratio(const Walker* walker_new, const Walker* walker_old, int i) const;
    virtual double get_val(const Walker* walker) const;
    virtual double get_lapl_sum(const Walker* walker) const;

};

#endif	/* JASTROW_H */

