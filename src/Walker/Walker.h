/* 
 * File:   Walker.h
 * Author: jorgehog
 *
 * Created on 30. mars 2012, 16:50
 */

#ifndef WALKER_H
#define	WALKER_H

using namespace arma;

class Walker {
protected:
    int n_p;
    int n2;
    int dim;

    bool is_murdered;

public:

    Walker(int n_p, int dim, bool do_init = true);

    double spatial_ratio;
    double value;
    double lapl_sum;
    double E;

    mat r;
    mat r_rel;

    mat qforce;

    mat spatial_grad;
    mat jast_grad;
    mat inv;

    mat phi;
    arma::field<mat> dell_phi;
    cube dJ;

    rowvec r2;



    void calc_r_i2(int i);

    void calc_r_i2() {
        for (int i = 0; i < n_p; i++) {
            this->calc_r_i2(i);
        }
    };

    double abs_relative(int i, int j) const;

    void make_rel_matrix();

    double get_r_i2(int i) const {
        return r2(i);
    };

    void kill() {
        is_murdered = true;
    }

    bool is_dead() {
        return is_murdered;
    }

    bool is_alive() {
        return !is_murdered;
    }

    void ressurect() {
        is_murdered = false;
    }

    void set_E(double E) {
        this->E = E;
    }

    double get_E() const {
        return E;
    }

    void print(std::string header = "----");


};

#endif	/* WALKER_H */
