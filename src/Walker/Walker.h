/* 
 * File:   Walker.h
 * Author: jorgehog
 *
 * Created on 30. mars 2012, 16:50
 */

#ifndef WALKER_H
#define	WALKER_H

class Walker {
protected:
    int n_p;
    int n2;
    int dim;
    
    bool is_murdered;

public:

    Walker();
    Walker(int n_p, int dim, bool do_init = true);
    
    double spatial_ratio;
    double lapl_sum;
    double E;

    arma::mat r;
    arma::mat r_rel;

    arma::mat qforce;

    arma::mat spatial_grad;
    arma::mat jast_grad;
    arma::mat inv;

    arma::mat phi;
    arma::field<arma::mat> dell_phi;
    arma::cube dJ;

    arma::rowvec r2;



    void calc_r_i2(int i);

    void calc_r_i2() {
        for (int i = 0; i < n_p; i++) {
            this->calc_r_i2(i);
        }
    }

    double calc_r_rel(int i, int j) const;

    void make_rel_matrix();
    
    void send_soul(int source);
    void recv_soul(int root);

    double get_r_i2(int i) const {
        return r2(i);
    }
    
    double get_r_i(int i) const {
        return sqrt(get_r_i2(i));
    }

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
