/* 
 * File:   Simple.h
 * Author: jorgmeister
 *
 * Created on October 12, 2012, 2:38 PM
 */

#ifndef SIMPLE_H
#define	SIMPLE_H

class Simple : public Diffusion {
public:
    Simple(int n_p, int dim, double timestep, long random_seed, double D = 0.5);

    virtual double get_new_pos(const Walker* walker, int i, int j) {
        return Diffusion::get_new_pos(walker, i, j);
    }

    virtual double get_g_ratio(const Walker* walker_post, const Walker* walker_pre) const {
        return 1;
    }

};

#endif	/* SIMPLE_H */

