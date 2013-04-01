/* 
 * File:   No_Jastrow.h
 * Author: jorgmeister
 *
 * Created on October 12, 2012, 3:23 PM
 */

#ifndef NO_JASTROW_H
#define	NO_JASTROW_H


/*!
 * \brief Class loaded when no correlation factor is used.
 */
class No_Jastrow : public Jastrow {
protected:

    double get_parameter(int n) {
        return 0;
    }
    
    void set_parameter(double param, int n) {};

    double get_variational_derivative(const Walker* walker, int n) {
        return 0;
    }

public:

    No_Jastrow();

    void get_grad(Walker* walker) const {};
    void get_grad(const Walker* walker_pre, Walker* walker_post, int i) const {};

    void initialize() {};
    
    void get_dJ_matrix(Walker* walker, int i) const {};

    double get_j_ratio(const Walker* walker_post, const Walker* walker_pre, int i) const {
        return 1;
    }

    double get_val(const Walker* walker) const {
        return 1;
    }

    double get_lapl_sum(Walker* walker) const {
        return 0;
    }

};

#endif	/* NO_JASTROW_H */
