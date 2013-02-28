/* 
 * File:   Diffusion.h
 * Author: jorgehog
 *
 * Created on 16. april 2012, 14:03
 */

#ifndef DIFFUSION_H
#define	DIFFUSION_H

class Diffusion {
protected:
    int n_p;
    int dim;

    QMC* qmc;

    double timestep;
    double D;
    long random_seed;

    double std;


public:

    Diffusion(int n_p, int dim, double timestep, seed_type random_seed, double D);

    virtual double get_new_pos(const Walker* walker, int i, int j);

    virtual double get_g_ratio(const Walker* walker_post, const Walker* walker_pre) const = 0;

    double get_GBfunc(double E_x, double E_y, double E_T) const {
        return exp(-(0.5 * (E_x + E_y) - E_T) * timestep);
    }

    double call_RNG();

    void set_qmc_ptr(QMC* qmc) {
        this->qmc = qmc;
    }

    void set_dt(double dt) {
        this->timestep = dt;
        this->std = sqrt(2 * D * dt);
    }

    double get_dt() const {
        return timestep;
    }

};


#endif	/* DIFFUSION_H */

