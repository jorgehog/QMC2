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
    
    double slater_ratio;
    double value;
    double lapl_sum;

    mat r;
    mat r_rel;

    mat qforce;

    mat spatial_grad;
    mat jast_grad;
    mat inv;

    rowvec r2;

    double get_r_i2(int i) const;

    void calc_r_i2(int i);
    void calc_r_i2();

    double abs_relative(int i, int j) const;

    void make_rel_matrix();

    bool is_singular() const;

    bool check_bad_qforce();

    void kill() {
        is_murdered = true;
    }

    bool is_dead() {
        return is_murdered;
    }

    void ressurect() {
        is_murdered = false;
    }


};

#endif	/* WALKER_H */
