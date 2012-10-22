/* 
 * File:   Pade_Jastrow.h
 * Author: jorgmeister
 *
 * Created on October 12, 2012, 2:38 PM
 */

#ifndef PADE_JASTROW_H
#define	PADE_JASTROW_H

class Pade_Jastrow : public Jastrow {
protected:
    double beta;
    arma::mat a;

    virtual double get_variational_derivative(const Walker* walker, int n) const;

    virtual void set_parameter(double param, int n) {
        beta = param;
    }

    virtual double get_parameter(int n) {
        return beta;
    }

public:

    Pade_Jastrow(GeneralParams &, VariationalParams &);

    virtual void initialize();

    virtual void get_grad(Walker* walker) const;
    virtual void get_grad(const Walker* walker_pre, Walker* walker_post, int i) const;
    virtual void get_dJ_matrix(Walker *walker, int i) const;


    virtual double get_j_ratio(const Walker* walker_new, const Walker* walker_old, int i) const;
    virtual double get_val(const Walker* walker) const;
    virtual double get_lapl_sum(const Walker* walker) const;

};

#endif	/* PADE_JASTROW_H */

