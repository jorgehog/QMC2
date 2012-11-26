/* 
 * File:   Sampling.h
 * Author: jorgehog
 *
 * Created on 15. juni 2012, 18:44
 */

#ifndef SAMPLING_H
#define	SAMPLING_H

class Sampling {
protected:
    int n_p;
    int n2;
    int dim;

    int start;
    int end;

    Diffusion* diffusion;
    QMC* qmc;

public:
    Sampling(int n_p, int dim);
    Sampling();

    void update_pos(const Walker* walker_pre, Walker* walker_post, int particle) const;
    virtual void update_necessities(const Walker* walker_pre, Walker* walker_post, int particle) const = 0;
    virtual void update_walker(Walker* walker_pre, const Walker* walker_post, int particle) const = 0;

    void set_trial_pos(Walker* walker, bool set_pos = true);
    void set_trial_states(Walker* walker);
    
    virtual void get_necessities(Walker* walker) = 0;
    virtual void calculate_energy_necessities(Walker* walker) const = 0;

    virtual void copy_walker(const Walker* parent, Walker* child) const = 0;
    virtual void reset_walker(const Walker* walker_pre, Walker* walker_post, int particle) const = 0;

    virtual double get_g_ratio(const Walker* walker_post, const Walker* walker_pre) const {
        return diffusion->get_g_ratio(walker_post, walker_pre);
    }

    double get_branching_Gfunc(double E_x, double E_y, double E_T) const {
        return diffusion->get_GBfunc(E_x, E_y, E_T);
    }

    double get_spatialjast_ratio(const Walker* walker_post, const Walker* walker_pre, int particle) const {
        return walker_post->spatial_ratio * qmc->get_jastrow_ptr()->get_j_ratio(walker_post, walker_pre, particle);
    }

    void set_qmc_ptr(QMC* qmc) {
        diffusion->set_qmc_ptr(qmc);
        this->qmc = qmc;
    }

    void set_dt(double dt) {
        diffusion->set_dt(dt);
    }

    double get_dt() const {
        return diffusion->get_dt();
    }

    double call_RNG() {
        return diffusion->call_RNG();
    }

    void set_spin_state(int start, int end) {
        this->start = start;
        this->end = end;
    }

};

#endif	/* SAMPLING_H */
