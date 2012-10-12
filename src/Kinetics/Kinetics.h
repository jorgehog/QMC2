/* 
 * File:   Kinetics.h
 * Author: jorgehog
 *
 * Created on 13. april 2012, 17:45
 */

#ifndef KINETICS_H
#define	KINETICS_H

class Kinetics {
protected:
    int n_p;
    int n2;
    int dim;

    QMC *qmc;

public:

    Kinetics();
    Kinetics(int n_p, int dim);

    virtual double get_KE(const Walker* walker) = 0;
    virtual void get_QF(Walker* walker) = 0;

    virtual void get_necessities_IS(Walker* walker) const = 0;
    //BF necessities are independent of kinetics

    virtual void update_walker_IS(Walker* walker_pre, const Walker* walker_post, int particle) const = 0;
    //NOTHING TO UPDATE FOR BF

    virtual double get_spatial_ratio_IS(const Walker* walker_post, const Walker* walker_pre, int particle) const = 0;
    //BF ratio is indep. of kinetics

    virtual void calculate_energy_necessities(Walker* walker) const = 0;

    virtual void update_necessities_IS(const Walker* walker_pre, Walker* walker_post, int particle) const = 0;
    //NO BF NECESSITIES

    virtual void copy_walker_IS(const Walker* parent, Walker* child) const = 0;
    virtual void copy_walker_BF(const Walker* parent, Walker* child) const = 0;

    virtual void reset_walker_IS(const Walker* walker_pre, Walker* walker_post, int particle) const = 0;
    //NOTHING TO RESET FOR BF

    void set_qmc_ptr(QMC* qmc) {
        this->qmc = qmc;
    }
};

#endif	/* KINETICS_H */

