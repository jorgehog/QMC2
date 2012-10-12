/* 
 * File:   Numerical.h
 * Author: jorgmeister
 *
 * Created on October 12, 2012, 2:39 PM
 */

#ifndef NUMERICAL_H
#define	NUMERICAL_H

class Numerical : public Kinetics {
protected:
    double h, h2;

    Walker* wfplus;
    Walker* wfminus;

public:
    Numerical();
    Numerical(GeneralParams &);

    virtual double get_KE(const Walker* walker);
    virtual void get_QF(Walker* walker);

    virtual void get_necessities_IS(Walker* walker) const;

    virtual void update_walker_IS(Walker* walker_pre, const Walker* walker_post, int particle) const;

    virtual double get_spatial_ratio_IS(const Walker* walker_post, const Walker* walker_pre, int particle) const;

    virtual void calculate_energy_necessities(Walker* walker) const;

    virtual void update_necessities_IS(const Walker* walker_pre, Walker* walker_post, int particle) const;

    virtual void copy_walker_IS(const Walker* parent, Walker* child) const;
    virtual void copy_walker_BF(const Walker* parent, Walker* child) const;

    virtual void reset_walker_IS(const Walker* walker_pre, Walker* walker_post, int particle) const;

};

#endif	/* NUMERICAL_H */

