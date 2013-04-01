/* 
 * File:   Pade_Jastrow.h
 * Author: jorgmeister
 *
 * Created on October 12, 2012, 2:38 PM
 */

#ifndef PADE_JASTROW_H
#define	PADE_JASTROW_H


/*!
 * \brief The Pade Jastrow factor with a single variational parameter.
 */
class Pade_Jastrow : public Jastrow {
protected:
    double beta; //!< The variational parameter
    arma::mat a; //!< The spin-dependent variables taking care of the cusp condition.

    double get_variational_derivative(const Walker* walker, int n);

    void set_parameter(double param, int n) {
        beta = param;
    }

    double get_parameter(int n) {
        return beta;
    }

public:

    Pade_Jastrow(GeneralParams &, VariationalParams &);

    /*!
     * In case of Pade Jastrow, initializing means seting up the a matrix.
     */
    void initialize();

    void get_grad(Walker* walker) const;
    void get_grad(const Walker* walker_pre, Walker* walker_post, int i) const;
    void get_dJ_matrix(Walker *walker, int i) const;


    double get_j_ratio(const Walker* walker_new, const Walker* walker_old, int i) const;
    double get_val(const Walker* walker) const;
    double get_lapl_sum(Walker* walker) const;

};

#endif	/* PADE_JASTROW_H */

