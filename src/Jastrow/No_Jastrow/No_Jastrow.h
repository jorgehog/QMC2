/* 
 * File:   No_Jastrow.h
 * Author: jorgmeister
 *
 * Created on October 12, 2012, 3:23 PM
 */

#ifndef NO_JASTROW_H
#define	NO_JASTROW_H

class No_Jastrow : public Jastrow {
protected:

    virtual double get_parameter(int n) {
        return 0;
    }
    
    virtual void set_parameter(double param, int n) {};

    virtual double get_variational_derivative(const Walker* walker, int n) const {
        return 0;
    }

public:

    No_Jastrow();

    virtual void get_grad(Walker* walker) const {};
    virtual void get_grad(const Walker* walker_pre, Walker* walker_post, int i) const {};

    virtual void initialize() {};
    
    virtual void get_dJ_matrix(Walker* walker, int i) const {};

    virtual double get_j_ratio(const Walker* walker_post, const Walker* walker_pre, int i) const {
        return 1;
    }

    virtual double get_val(const Walker* walker) const {
        return 1;
    }

    virtual double get_lapl_sum(const Walker* walker) const {
        return 0;
    }

};

#endif	/* NO_JASTROW_H */
